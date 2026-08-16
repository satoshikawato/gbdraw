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
const visibilityActionsSource = await readFile(
  join(sourceDir, 'feature-editor', 'visibility-actions.js'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'feature-editor', 'visibility-actions.js'),
  visibilityActionsSource,
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
assert.doesNotMatch(visibilityActionsSource, /flushActiveResult|serializeCleanSvg|results\.value\[[^\]]+\]\s*=/);

const ref = (value) => ({ value });
const featureA = { svg_id: 'feature-a', type: 'CDS', label: 'A' };
const featureB = { svg_id: 'feature-b', type: 'CDS', label: 'B' };
const extractedFeatures = ref([featureA, featureB]);
const orthogroups = ref([]);
const featureVisibilityManualRules = [];
const featureVisibilityOverrides = {};
const clickedFeature = ref({ svg_id: 'feature-a', featureVisibility: 'default' });
const featureVisibilityScopeDialog = {};
const selectedResultIndex = ref(0);
const resultGenerationKey = ref('generation-1');
const appliedPreviewChanges = [];
const labelLayoutDirtyReason = ref('');

const actions = createFeatureVisibilityActions({
  state: {
    clickedFeature,
    extractedFeatures,
    orthogroups,
    featureVisibilityManualRules,
    featureVisibilityRules: ref([]),
    featureVisibilityOverrides,
    featureVisibilitySelectorCache: {},
    featureVisibilityScopeDialog,
    labelLayoutDirtyReason,
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

assert.equal(await command.revert(), true);
assert.deepEqual(featureVisibilityOverrides, {});
assert.equal(appliedPreviewChanges.length, 2);
assert.deepEqual(
  appliedPreviewChanges[1].changes.map((change) => [change.featureId, change.mode]),
  [['feature-a', 'on'], ['feature-b', 'on']]
);

assert.equal(actions.setFeatureVisibility(featureA, 'off', {
  triggerReflow: false,
  scope: { id: 'feature' }
}), true);
assert.equal(featureVisibilityOverrides['feature-a'], 'off');
assert.equal(appliedPreviewChanges.length, 3);
assert.deepEqual(appliedPreviewChanges.at(-1).changes, [{ featureId: 'feature-a', mode: 'off' }]);
delete featureVisibilityOverrides['feature-a'];

const productFeatureA = {
  svg_id: 'product-feature-a',
  type: 'CDS',
  product: 'shared product',
  qualifiers: { protein_id: ['WP_PRODUCT_A.1'] }
};
const productFeatureB = {
  svg_id: 'product-feature-b',
  type: 'CDS',
  product: 'shared product',
  qualifiers: { protein_id: ['WP_PRODUCT_B.1'] }
};
extractedFeatures.value = [productFeatureA, productFeatureB];
const productScope = {
  id: 'product',
  featureType: 'CDS',
  qualifier: 'product',
  value: 'shared product',
  label: 'Exact product: shared product'
};
assert.equal(actions.setFeatureVisibility(productFeatureA, 'off', {
  triggerReflow: false,
  scope: productScope
}), true);
assert.deepEqual(
  appliedPreviewChanges.at(-1).changes,
  [
    { featureId: 'product-feature-a', mode: 'off' },
    { featureId: 'product-feature-b', mode: 'off' }
  ]
);
assert.equal(featureVisibilityManualRules.at(-1).qualifier, 'product');
actions.setFeatureVisibility(productFeatureA, 'default', {
  triggerReflow: false,
  scope: productScope
});
assert.deepEqual(
  appliedPreviewChanges.at(-1).changes,
  [
    { featureId: 'product-feature-a', mode: 'on' },
    { featureId: 'product-feature-b', mode: 'on' }
  ]
);
assert.equal(featureVisibilityManualRules.length, 0);

const proteinScope = {
  id: 'protein_id',
  featureType: 'CDS',
  qualifier: 'protein_id',
  value: 'WP_PRODUCT_B.1',
  label: 'Exact protein ID: WP_PRODUCT_B.1'
};
assert.equal(actions.setFeatureVisibility(productFeatureB, 'off', {
  triggerReflow: false,
  scope: proteinScope
}), true);
assert.deepEqual(
  appliedPreviewChanges.at(-1).changes,
  [{ featureId: 'product-feature-b', mode: 'off' }]
);
assert.equal(featureVisibilityManualRules.at(-1).qualifier, 'protein_id');
featureVisibilityManualRules.splice(0);

extractedFeatures.value = [productFeatureA, productFeatureB];
assert.equal(actions.addFeatureVisibilityRule(), true);
assert.equal(labelLayoutDirtyReason.value, 'feature-visibility-rule-add');
assert.deepEqual(
  appliedPreviewChanges.at(-1).changes.map((change) => change.mode),
  ['on', 'on']
);
assert.equal(actions.setFeatureVisibilityRuleField(0, 'featureType', 'CDS'), true);
assert.equal(actions.setFeatureVisibilityRuleField(0, 'value', '^shared product$'), true);
assert.deepEqual(
  appliedPreviewChanges.at(-1).changes.map((change) => change.mode),
  ['off', 'off']
);
assert.equal(actions.setFeatureVisibilityRuleField(0, 'action', 'show'), true);
assert.deepEqual(
  appliedPreviewChanges.at(-1).changes.map((change) => change.mode),
  ['on', 'on']
);
assert.equal(actions.removeFeatureVisibilityRule(0), true);
assert.deepEqual(
  appliedPreviewChanges.at(-1).changes.map((change) => change.mode),
  ['on', 'on']
);

featureVisibilityManualRules.push(
  {
    id: 'ordered-off', source: 'manual', recordId: '*', featureType: 'CDS',
    qualifier: 'product', value: '^shared product$', action: 'off'
  },
  {
    id: 'ordered-on', source: 'manual', recordId: '*', featureType: 'CDS',
    qualifier: 'product', value: '^shared product$', action: 'show'
  }
);
actions.reconcileFeatureVisibility();
assert.deepEqual(
  appliedPreviewChanges.at(-1).changes.map((change) => change.mode),
  ['off', 'off']
);
assert.equal(actions.moveFeatureVisibilityRuleDown(0), true);
assert.deepEqual(
  appliedPreviewChanges.at(-1).changes.map((change) => change.mode),
  ['on', 'on']
);
featureVisibilityManualRules.splice(0);
extractedFeatures.value = [featureA, featureB];

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

const previewCallCountBeforeStaleCommand = appliedPreviewChanges.length;
resultGenerationKey.value = 'generation-2';
assert.equal(await command.apply(), false);
assert.deepEqual(featureVisibilityOverrides, {});
assert.equal(appliedPreviewChanges.length, previewCallCountBeforeStaleCommand);

console.log('feature visibility action tests passed');
