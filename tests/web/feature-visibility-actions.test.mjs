import assert from 'node:assert/strict';
import { createFeatureVisibilityActions } from '../../gbdraw/web/js/app/feature-editor/visibility-actions.js';
import {
  requestFeatureVisibilityRules,
  setFeatureVisibilityOverride
} from '../../gbdraw/web/js/app/feature-visibility.js';
import { createRulePreparation } from '../../gbdraw/web/js/app/rule-matching.js';
import { evaluatePythonRules } from './helpers/python-rule-evaluator.mjs';

const ref = (value) => ({ value });
// The live projection reads Python's rule matches (R4) and the committed
// request's feature types.
const rulePreparationFor = (state) => createRulePreparation({
  state,
  evaluate: evaluatePythonRules,
  visibilityRules: () => requestFeatureVisibilityRules(state.featureVisibilityManualRules)
});
const committedRequest = (selectedFeaturesSet) => () => ({ diagramOptions: { selectedFeaturesSet } });
// Per-feature visibility is the identity row of the feature (design Q4).
const identity = (id) => ({ scope: 'circular', record_key: 'record-1', biological_feature_id: `bio-${id}` });
const modes = (overrides) => Object.fromEntries(Object.values(overrides)
  .filter((row) => row.featureVisibility !== null)
  .map((row) => [row.biologicalFeatureId.replace(/^bio-/, ''), row.featureVisibility]));
const clear = (overrides) => Object.keys(overrides).forEach((key) => delete overrides[key]);
const featureA = { svg_id: 'feature-a', type: 'CDS', label: 'A', ...identity('feature-a') };
const featureB = { svg_id: 'feature-b', type: 'CDS', label: 'B', ...identity('feature-b') };
const extractedFeatures = ref([featureA, featureB]);
const orthogroups = ref([]);
const featureVisibilityOverrides = {};
const clickedFeature = ref({ svg_id: 'feature-a', featureVisibility: 'default' });
const featureVisibilityScopeDialog = {};
const selectedResultIndex = ref(0);
const resultGenerationKey = ref('generation-1');
const appliedPreviewChanges = [];
// F-3: each visibility edit hands the feature's label to the label owner.
const labelVisibilityCalls = [];

const actionState = {
  clickedFeature,
  extractedFeatures,
  orthogroups,
  featureVisibilityManualRules: [],
  featureVisibilityRules: ref([]),
  featureOverrides: featureVisibilityOverrides,
  featureVisibilityScopeDialog,
  resultGenerationKey,
  results: ref([{ name: 'one.svg', content: '<svg></svg>' }]),
  selectedResultIndex,
  svgContainer: ref({
    querySelector: (selector) => (selector === 'svg' ? {} : null)
  })
};
const actions = createFeatureVisibilityActions({
  state: actionState,
  rulePreparation: rulePreparationFor(actionState),
  getCommittedRequest: committedRequest(['CDS']),
  featureSvgActions: {
    applyVisibilityPreviewChanges: (changes, options = {}) => {
      appliedPreviewChanges.push({ changes, reason: options.reason });
      return true;
    }
  },
  labelActions: {
    applyFeatureVisibilityToLabels: (options = {}) => {
      labelVisibilityCalls.push(options.reflow !== false);
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
assert.deepEqual(modes(featureVisibilityOverrides), {
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
assert.deepEqual(labelVisibilityCalls, [
  true,
  true
]);
assert.deepEqual(
  appliedPreviewChanges[1].changes.map((change) => [change.featureId, change.mode]),
  [['feature-a', 'on'], ['feature-b', 'on']]
);

assert.equal(actions.setFeatureVisibility(featureA, 'off', {
  triggerReflow: false,
  scope: { id: 'feature' }
}), true);
assert.equal(modes(featureVisibilityOverrides)['feature-a'], 'off');
assert.equal(appliedPreviewChanges.length, 3);
assert.deepEqual(labelVisibilityCalls.at(-1), false,
  'the label follows the feature even when the caller declines the reflow');
assert.deepEqual(
  appliedPreviewChanges[2].changes.map((change) => [change.featureId, change.mode]),
  [['feature-a', 'off']]
);
clear(featureVisibilityOverrides);

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
await actions.updateClickedFeatureVisibility('off');
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
await actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, false);
clear(featureVisibilityOverrides);

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
await actions.updateClickedFeatureVisibility('off');
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
await actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, false);
clear(featureVisibilityOverrides);

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
await actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, false);
clear(featureVisibilityOverrides);

extractedFeatures.value = [strictTrigger, duplicateA];
orthogroups.value = [orthogroups.value[0], { ...orthogroups.value[0] }];
featureVisibilityScopeDialog.show = false;
featureVisibilityScopeDialog.scopes = [];
await actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, false);
clear(featureVisibilityOverrides);

const previewChangeCountBeforeStaleApply = appliedPreviewChanges.length;
resultGenerationKey.value = 'generation-2';
assert.equal(await command.apply(), false);
assert.deepEqual(featureVisibilityOverrides, {});
assert.equal(appliedPreviewChanges.length, previewChangeCountBeforeStaleApply);

// FE-04: an Exact product hide is committed as an editor qualifier rule. The
// reconcile that History runs after Undo/Redo resolves it with the action's
// resolver (Python's matches), so the feature stays hidden. R-2: every rule,
// also one typed in the Features panel or loaded from a TSV, hides live as
// Generate hides.
{
  const nd1 = { svg_id: 'nd1', type: 'CDS', qualifiers: { product: ['NADH dehydrogenase subunit 1'] }, ...identity('nd1') };
  const nd1Case = {
    svg_id: 'nd1-case', type: 'CDS', qualifiers: { product: ['nadh DEHYDROGENASE subunit 1'] }, ...identity('nd1-case')
  };
  const nd2 = { svg_id: 'nd2', type: 'CDS', qualifiers: { product: ['NADH dehydrogenase subunit 2'] }, ...identity('nd2') };
  const geneRna = {
    svg_id: 'nd1-rna', type: 'tRNA', qualifiers: { product: ['NADH dehydrogenase subunit 1'] }, ...identity('nd1-rna')
  };
  const manualRules = [];
  const overrides = {};
  const reconciled = [];
  const scopeDialog = {};
  const productState = {
    clickedFeature: ref({ svg_id: nd1.svg_id, featureVisibility: 'default', feat: nd1 }),
    extractedFeatures: ref([nd1, nd1Case, nd2, geneRna]),
    orthogroups: ref([]),
    featureVisibilityManualRules: manualRules,
    featureVisibilityRules: ref([]),
    featureOverrides: overrides,
    featureVisibilityScopeDialog: scopeDialog,
    resultGenerationKey: ref('generation-1'),
    results: ref([{ name: 'one.svg', content: '<svg></svg>' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: (selector) => (selector === 'svg' ? {} : null) })
  };
  const productPreparation = rulePreparationFor(productState);
  const productActions = createFeatureVisibilityActions({
    state: productState,
    rulePreparation: productPreparation,
    getCommittedRequest: committedRequest(['CDS', 'tRNA']),
    featureSvgActions: {
      applyVisibilityPreviewChanges: (changes) => {
        reconciled.push(Object.fromEntries(changes.map((change) => [change.featureId, change.mode])));
        return true;
      }
    },
    previewRuntime: { selectResult: () => true }
  });
  await productActions.updateClickedFeatureVisibility('off');
  assert.equal(scopeDialog.show, true);
  assert.equal(await productActions.handleFeatureVisibilityScopeChoice('product'), true);
  const hidden = { nd1: 'off', 'nd1-case': 'off', nd2: 'on', 'nd1-rna': 'on' };
  assert.deepEqual(reconciled.at(-1), hidden);
  assert.equal(manualRules.length, 1);
  assert.equal(productActions.reconcileFeatureVisibility(), true);
  assert.deepEqual(reconciled.at(-1), hidden);
  setFeatureVisibilityOverride(overrides, nd1Case, 'on');
  productActions.reconcileFeatureVisibility();
  assert.equal(reconciled.at(-1)['nd1-case'], 'on', 'a per-feature override takes precedence');
  clear(overrides);
  manualRules.splice(0, manualRules.length, {
    ...manualRules[0], source: 'file', recordId: '*', featureType: '*', qualifier: 'Product', value: 'subunit 2$'
  });
  assert.equal(await productPreparation.prepareDrawn(), true);
  productActions.reconcileFeatureVisibility();
  assert.deepEqual(reconciled.at(-1), { nd1: 'on', 'nd1-case': 'on', nd2: 'off', 'nd1-rna': 'on' },
    'a loaded rule hides live as Generate does');
}

console.log('feature visibility action tests passed');
