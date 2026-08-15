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
await writeFile(join(tempDir, 'app', 'feature-utils.js'), await readFile(join(sourceDir, 'feature-utils.js'), 'utf8'), 'utf8');
await writeFile(
  join(tempDir, 'services', 'text-download.js'),
  await readFile(join(sourceDir, '..', 'services', 'text-download.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'feature-identity.js'),
  await readFile(join(sourceDir, '..', 'services', 'feature-identity.js'), 'utf8'),
  'utf8'
);

const { createFeatureVisibilityActions } = await import(
  pathToFileURL(join(tempDir, 'app', 'feature-editor', 'visibility-actions.js'))
);

const ref = (value) => ({ value });
const featureA = { svg_id: 'feature-a', type: 'CDS', label: 'A' };
const featureB = { svg_id: 'feature-b', type: 'CDS', label: 'B' };
const extractedFeatures = ref([featureA, featureB]);
const orthogroups = ref([]);
const featureVisibilityOverrides = {};
const clickedFeature = ref({ svg_id: 'feature-a', featureVisibility: 'default' });
const featureVisibilityScopeDialog = {};
const selectedResultIndex = ref(0);
const resultGenerationKey = ref('generation-1');
const appliedPreviewChanges = [];
const flushes = [];

const actions = createFeatureVisibilityActions({
  state: {
    clickedFeature,
    extractedFeatures,
    orthogroups,
    featureVisibilityManualRules: [],
    featureVisibilityRules: ref([]),
    featureVisibilityOverrides,
    featureVisibilitySelectorCache: {},
    featureVisibilityScopeDialog,
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

assert.equal(actions.setFeatureVisibility(featureA, 'off', {
  triggerReflow: false,
  scope: { id: 'feature' }
}), true);
assert.equal(featureVisibilityOverrides['feature-a'], 'off');
assert.deepEqual(flushes, [{ force: true }, { force: true }, {}]);
delete featureVisibilityOverrides['feature-a'];

const sourceIdFeature = {
  svg_id: 'feature-source-id',
  type: 'CDS',
  proteinId: 'h_aaaaaaaaaaaaaaaaaaaaaaaaaa',
  sourceProteinId: 'WP_012345678.1'
};
clickedFeature.value = {
  svg_id: sourceIdFeature.svg_id,
  featureVisibility: 'default',
  feat: sourceIdFeature
};
actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, true);
assert.ok(featureVisibilityScopeDialog.scopes.some((scope) => (
  scope.id === 'protein_id' &&
  scope.value === 'WP_012345678.1' &&
  scope.label === 'Exact protein ID: WP_012345678.1'
)));
assert.ok(featureVisibilityScopeDialog.scopes.every((scope) => !scope.label.includes('h_')));

const runtimeOnlyFeature = {
  svg_id: 'feature-runtime-only',
  type: 'CDS',
  proteinId: 'h_bbbbbbbbbbbbbbbbbbbbbbbbbb'
};
clickedFeature.value = {
  svg_id: runtimeOnlyFeature.svg_id,
  featureVisibility: 'default',
  feat: runtimeOnlyFeature
};
featureVisibilityScopeDialog.show = false;
actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, false);
delete featureVisibilityOverrides[runtimeOnlyFeature.svg_id];

const reversedFeature = {
  svg_id: 'display-feature_record_3',
  stable_feature_id: 'source-feature',
  record_idx: 2,
  type: 'CDS',
  orthogroupId: 'og_reverse'
};
const groupedFeature = {
  svg_id: 'grouped-feature_record_1',
  stable_feature_id: 'grouped-source-feature',
  record_idx: 0,
  type: 'CDS',
  orthogroupId: 'og_reverse'
};
extractedFeatures.value.push(reversedFeature, groupedFeature);
orthogroups.value = [{
  id: 'og_reverse',
  members: [
    {
      recordIndex: 2,
      featureSvgId: 'source-feature',
      stableFeatureSvgId: 'source-feature',
      renderedFeatureSvgId: 'display-feature_record_3'
    },
    {
      recordIndex: 0,
      featureSvgId: 'grouped-source-feature',
      stableFeatureSvgId: 'grouped-source-feature',
      renderedFeatureSvgId: 'grouped-feature_record_1'
    }
  ]
}];
clickedFeature.value = {
  svg_id: reversedFeature.svg_id,
  featureVisibility: 'default',
  feat: reversedFeature
};
actions.updateClickedFeatureVisibility('off');
const reversedGroupScope = featureVisibilityScopeDialog.scopes.find((scope) => scope.id === 'orthogroup');
assert.ok(reversedGroupScope);
assert.deepEqual(
  reversedGroupScope.features.map((feature) => feature.svg_id),
  ['display-feature_record_3', 'grouped-feature_record_1']
);

const strictTrigger = {
  svg_id: 'strict-trigger-rendered',
  stable_feature_id: 'strict-trigger-source',
  record_idx: 2,
  type: 'CDS',
  orthogroupId: 'og_strict'
};
const wrongRecordFeature = {
  svg_id: 'shared-rendered',
  stable_feature_id: 'shared-source',
  record_idx: 0,
  type: 'CDS',
  orthogroupId: 'og_strict'
};
extractedFeatures.value = [strictTrigger, wrongRecordFeature];
orthogroups.value = [{
  id: 'og_strict',
  members: [
    {
      recordIndex: 2,
      featureSvgId: 'strict-trigger-source',
      stableFeatureSvgId: 'strict-trigger-source',
      renderedFeatureSvgId: 'strict-trigger-rendered'
    },
    {
      recordIndex: 1,
      featureSvgId: 'shared-source',
      stableFeatureSvgId: 'shared-source',
      renderedFeatureSvgId: 'shared-rendered'
    }
  ]
}];
clickedFeature.value = {
  svg_id: strictTrigger.svg_id,
  featureVisibility: 'default',
  feat: strictTrigger
};
featureVisibilityScopeDialog.show = false;
featureVisibilityScopeDialog.scopes = [];
actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, false);
delete featureVisibilityOverrides[strictTrigger.svg_id];

const duplicateA = {
  svg_id: 'duplicate-rendered',
  stable_feature_id: 'duplicate-source',
  record_idx: 1,
  type: 'CDS',
  orthogroupId: 'og_strict'
};
const duplicateB = { ...duplicateA, label: 'duplicate metadata row' };
extractedFeatures.value = [strictTrigger, duplicateA, duplicateB];
orthogroups.value[0].members[1] = {
  recordIndex: 1,
  featureSvgId: 'duplicate-source',
  stableFeatureSvgId: 'duplicate-source',
  renderedFeatureSvgId: 'duplicate-rendered'
};
featureVisibilityScopeDialog.show = false;
featureVisibilityScopeDialog.scopes = [];
actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, false);
delete featureVisibilityOverrides[strictTrigger.svg_id];

extractedFeatures.value = [strictTrigger, duplicateA];
orthogroups.value = [orthogroups.value[0], { ...orthogroups.value[0] }];
featureVisibilityScopeDialog.show = false;
featureVisibilityScopeDialog.scopes = [];
actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, false);
delete featureVisibilityOverrides[strictTrigger.svg_id];

flushes.length = 0;
resultGenerationKey.value = 'generation-2';
assert.equal(await command.apply(), false);
assert.deepEqual(featureVisibilityOverrides, {});
assert.equal(appliedPreviewChanges.length, 2);
assert.deepEqual(flushes, []);

console.log('feature visibility action tests passed');
