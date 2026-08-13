import assert from 'node:assert/strict';
import test from 'node:test';

import {
  FEATURE_LABEL_SCOPE_ORDER,
  featureLabelStateFingerprint,
  planFeatureLabelChange,
  resolveFeatureLabelPlan,
  stableFeatureTargetKey
} from '../../gbdraw/web/js/app/feature-editor/label-scope-plan.js';
import { resolveFeatureLabelViewModel } from '../../gbdraw/web/js/app/feature-editor/label-view-model.js';

const hashes = {
  a: 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa',
  b: 'fi1_bbbbbbbbbbbbbbbbbbbbbbbbbb',
  c: 'fi1_cccccccccccccccccccccccccc'
};

const biological = (recordKey, id, instanceHash, product) => ({
  recordKey,
  biologicalFeatureId: id,
  instanceHash,
  type: 'CDS',
  product,
  qualifiers: { product: [product] }
});

const catalog = {
  schema: 3,
  items: [{
    resultKey: 'result-a', resultIndex: 0, resultName: 'a.svg', recordKeys: ['record-a'],
    biologicalFeatures: [biological('record-a', 'feature-a', hashes.a, 'Kinase')],
    features: [{ recordKey: 'record-a', biologicalFeatureId: 'feature-a', svgId: 'svg-a' }]
  }, {
    resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg', recordKeys: ['record-b'],
    biologicalFeatures: [
      biological('record-b', 'feature-b', hashes.b, 'Kinase'),
      biological('record-b', 'feature-c', hashes.c, 'Transporter')
    ],
    features: [
      { recordKey: 'record-b', biologicalFeatureId: 'feature-b', svgId: 'svg-b' },
      { recordKey: 'record-b', biologicalFeatureId: 'feature-c', svgId: 'svg-c' }
    ]
  }]
};

const feature = (resultIndex, featureIndex) => {
  const item = catalog.items[resultIndex];
  const source = item.biologicalFeatures[featureIndex];
  return {
    ...source,
    resultKey: item.resultKey,
    resultIndex,
    resultName: item.resultName,
    renderedFeatureSvgId: item.features[featureIndex].svgId
  };
};

const targetA = stableFeatureTargetKey(feature(0, 0));
const targetC = stableFeatureTargetKey(feature(1, 1));

const deepFreeze = (value) => {
  if (!value || typeof value !== 'object' || Object.isFrozen(value)) return value;
  Object.values(value).forEach(deepFreeze);
  return Object.freeze(value);
};

const baseSnapshot = () => ({
  catalog,
  presentationLabelsBySvgId: {
    'svg-a': { text: 'Shared rendered', sourceText: 'Kinase', labelKey: 'label-a' },
    'svg-b': { text: 'Shared rendered', sourceText: 'Kinase', labelKey: 'label-b' },
    'svg-c': { text: 'Shared rendered', sourceText: 'Transporter', labelKey: 'label-c' }
  },
  labelTextFeatureOverrides: {},
  labelTextBulkOverrides: {},
  labelTextFeatureOverrideSources: {},
  selectedFeatureKeys: [targetA, targetC],
  exactScopeAvailable: true,
  resultGenerationKey: 7,
  documentEpoch: 3,
  styleRevision: 4,
  styleFingerprint: 'sf1_example'
});

const intent = (overrides = {}) => ({
  targetFeatureKey: targetA,
  newText: 'Renamed',
  source: 'popup',
  originResultKey: 'result-a',
  originResultIndex: 0,
  resultGenerationKey: 7,
  documentEpoch: 3,
  styleRevision: 4,
  styleFingerprint: 'sf1_example',
  ...overrides
});

test('label view model separates rendered text, source text, and persisted origins', () => {
  const inherited = resolveFeatureLabelViewModel({
    feature: biological('record-a', 'feature-a', hashes.a, 'Kinase'),
    renderedSvgIds: ['svg-a']
  });
  assert.equal(inherited.effectiveText, 'Kinase');
  assert.equal(inherited.origin, 'source-annotation');
  assert.equal(inherited.canReset, false);

  const grouped = resolveFeatureLabelViewModel({
    feature: biological('record-a', 'feature-a', hashes.a, 'Kinase'),
    renderedSvgIds: ['svg-a'],
    bulkOverrides: { Kinase: 'All kinases' }
  });
  assert.equal(grouped.effectiveText, 'All kinases');
  assert.equal(grouped.origin, 'source-group');

  const exact = resolveFeatureLabelViewModel({
    feature: biological('record-a', 'feature-a', hashes.a, 'Kinase'),
    renderedSvgIds: ['svg-a'],
    featureOverrides: { 'svg-a': 'Only this one' },
    featureOverrideSources: { 'svg-a': 'Kinase' }
  });
  assert.equal(exact.effectiveText, 'Only this one');
  assert.equal(exact.explicitValue, 'Only this one');
  assert.equal(exact.canReset, true);
});

test('planner returns deterministic rendered/source/selected/single scopes across Results', () => {
  const snapshot = baseSnapshot();
  const before = structuredClone(snapshot);
  deepFreeze(snapshot);
  const plan = planFeatureLabelChange(snapshot, deepFreeze(intent()));
  assert.equal(plan.status, 'needs-scope');
  assert.deepEqual(plan.candidates.map((candidate) => candidate.semanticScope), FEATURE_LABEL_SCOPE_ORDER);
  assert.deepEqual(plan.candidates.map((candidate) => candidate.targetCount), [3, 2, 2, 1]);
  assert.deepEqual(plan.candidates.map((candidate) => candidate.affectedResultCount), [2, 2, 2, 1]);
  assert.equal(plan.candidates[0].resultExtent, 'all-affected-results');
  assert.equal(plan.candidates.at(-1).resultExtent, 'current-result');
  assert.match(plan.candidates[0].label, /3 features in 2 outputs/);
  assert.equal(Object.isFrozen(plan), true);
  assert.deepEqual(snapshot, before, 'planning must not mutate its inputs');

  const source = plan.candidates.find((candidate) => (
    candidate.semanticScope === 'source-annotation-label'
  ));
  const resolved = resolveFeatureLabelPlan(plan, source.id);
  assert.equal(resolved.status, 'ready');
  assert.equal(resolved.semanticScope, 'source-annotation-label');
  assert.deepEqual(resolved.affectedResultKeys, ['result-a', 'result-b']);
  assert.deepEqual(resolved.targetsByResult.map((group) => group.featureKeys.length), [1, 1]);
  assert.equal(resolved.semanticAfter.durableLabelIntent.sourceText, 'Kinase');
});

test('explicit selection is ready while unsupported or stale requests fail closed', () => {
  const selected = planFeatureLabelChange(baseSnapshot(), intent({
    source: 'selection-toolbar',
    semanticScope: 'selected-features',
    targetFeatureKeys: [targetA, targetC]
  }));
  assert.equal(selected.status, 'ready');
  assert.equal(resolveFeatureLabelPlan(selected).semanticScope, 'selected-features');

  const legacy = baseSnapshot();
  legacy.exactScopeAvailable = false;
  const legacyPlan = planFeatureLabelChange(legacy, intent({ semanticScope: 'single' }));
  assert.equal(legacyPlan.status, 'invalid');
  assert.match(legacyPlan.diagnostics.map((entry) => entry.message).join(' '), /Generate to enable/);

  const stale = planFeatureLabelChange(baseSnapshot(), intent({ styleRevision: 5 }));
  assert.equal(stale.status, 'invalid');
  assert.equal(stale.diagnostics[0].code, 'stale-revision');

  const unsupported = planFeatureLabelChange(baseSnapshot(), intent({ semanticScope: 'matching-rule' }));
  assert.equal(unsupported.status, 'invalid');
  assert.equal(unsupported.diagnostics[0].code, 'unsupported-scope');
});

test('label-state fingerprint is stable by object key and distinguishes each owner', () => {
  const a = featureLabelStateFingerprint({
    labelTextFeatureOverrides: { b: 'B', a: 'A' },
    labelTextBulkOverrides: { Kinase: 'K' },
    labelTextFeatureOverrideSources: { a: 'Source' }
  });
  const b = featureLabelStateFingerprint({
    labelTextFeatureOverrides: { a: 'A', b: 'B' },
    labelTextBulkOverrides: { Kinase: 'K' },
    labelTextFeatureOverrideSources: { a: 'Source' }
  });
  assert.equal(a, b);
  assert.notEqual(a, featureLabelStateFingerprint({
    labelTextFeatureOverrides: { a: 'Changed', b: 'B' },
    labelTextBulkOverrides: { Kinase: 'K' },
    labelTextFeatureOverrideSources: { a: 'Source' }
  }));
});
