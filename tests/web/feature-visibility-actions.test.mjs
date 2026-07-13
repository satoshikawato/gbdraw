import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourceDir = join(repoRoot, 'gbdraw', 'web', 'js', 'app');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-feature-visibility-actions-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app', 'feature-editor'), { recursive: true });
await mkdir(join(tempDir, 'services'), { recursive: true });
await writeFile(
  join(tempDir, 'app', 'feature-editor', 'visibility-actions.js'),
  await readFile(join(sourceDir, 'feature-editor', 'visibility-actions.js'), 'utf8'),
  'utf8'
);
await writeFile(join(tempDir, 'app', 'feature-visibility.js'), await readFile(join(sourceDir, 'feature-visibility.js'), 'utf8'), 'utf8');
await writeFile(join(tempDir, 'app', 'feature-selector.js'), await readFile(join(sourceDir, 'feature-selector.js'), 'utf8'), 'utf8');
await writeFile(
  join(tempDir, 'services', 'text-download.js'),
  await readFile(join(sourceDir, '..', 'services', 'text-download.js'), 'utf8'),
  'utf8'
);

const { createFeatureVisibilityActions } = await import(
  pathToFileURL(join(tempDir, 'app', 'feature-editor', 'visibility-actions.js'))
);

const ref = (value) => ({ value });
const featureA = { svg_id: 'feature-a', type: 'CDS', label: 'A' };
const featureB = { svg_id: 'feature-b', type: 'CDS', label: 'B' };
const featureVisibilityOverrides = {};
const selectedResultIndex = ref(0);
const resultGenerationKey = ref('generation-1');
const appliedPreviewChanges = [];
const flushes = [];

const actions = createFeatureVisibilityActions({
  state: {
    clickedFeature: ref({ svg_id: 'feature-a', featureVisibility: 'default' }),
    extractedFeatures: ref([featureA, featureB]),
    orthogroups: ref([]),
    featureVisibilityManualRules: [],
    featureVisibilityRules: ref([]),
    featureVisibilityOverrides,
    featureVisibilitySelectorCache: {},
    featureVisibilityScopeDialog: {},
    labelLayoutDirtyReason: ref(''),
    resultGenerationKey,
    results: ref([{ name: 'one.svg', content: '<svg></svg>' }]),
    selectedResultIndex,
    svgContainer: ref({
      querySelector: (selector) => (selector === 'svg' ? {} : null)
    })
  },
  featureSvgActions: {
    applyVisibilityPreviewBySvgId: () => true,
    applyVisibilityPreviewChanges: (changes, options = {}) => {
      appliedPreviewChanges.push({ changes, reason: options.reason });
      return true;
    }
  },
  previewRuntime: {
    selectResult: (index) => {
      selectedResultIndex.value = index;
      return true;
    },
    flushActiveResult: (options = {}) => {
      flushes.push(options);
      return true;
    }
  }
});

const command = actions.buildSelectedFeaturesVisibilityCommand([featureA, featureB], 'off');
assert.ok(command);
assert.equal(await command.apply(), true);
assert.deepEqual(featureVisibilityOverrides, {
  'feature-a': 'off',
  'feature-b': 'off'
});
assert.equal(appliedPreviewChanges.length, 1);
assert.deepEqual(
  appliedPreviewChanges[0].changes.map((change) => [change.featureId, change.mode]),
  [['feature-a', 'off'], ['feature-b', 'off']]
);
assert.deepEqual(flushes, [{ force: true }]);

assert.equal(await command.revert(), true);
assert.deepEqual(featureVisibilityOverrides, {});
assert.equal(appliedPreviewChanges.length, 2);
assert.deepEqual(
  appliedPreviewChanges[1].changes.map((change) => [change.featureId, change.mode]),
  [['feature-a', 'on'], ['feature-b', 'on']]
);
assert.deepEqual(flushes, [{ force: true }, { force: true }]);

resultGenerationKey.value = 'generation-2';
assert.equal(await command.apply(), false);
assert.deepEqual(featureVisibilityOverrides, {});
assert.equal(appliedPreviewChanges.length, 2);
assert.deepEqual(flushes, [{ force: true }, { force: true }]);

console.log('feature visibility action tests passed');
