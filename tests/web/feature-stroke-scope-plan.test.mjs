import assert from 'node:assert/strict';
import test from 'node:test';

import {
  FEATURE_STROKE_SCOPE_ORDER,
  planFeatureStrokeChange,
  resolveFeatureStrokePlan,
  stableFeatureTargetKey
} from '../../gbdraw/web/js/app/feature-editor/stroke-scope-plan.js';
import { resolveFeatureStrokeViewModel } from '../../gbdraw/web/js/app/feature-editor/stroke-view-model.js';

const hashes = {
  a: 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa',
  b: 'fi1_bbbbbbbbbbbbbbbbbbbbbbbbbb',
  c: 'fi1_cccccccccccccccccccccccccc'
};

const feature = ({
  resultKey,
  resultIndex,
  resultName,
  recordKey,
  biologicalFeatureId,
  instanceHash,
  legendCaption = 'Core',
  svgId
}) => ({
  resultKey,
  resultIndex,
  resultName,
  recordKey,
  biologicalFeatureId,
  instanceHash,
  type: 'CDS',
  legendCaption,
  renderedFeatureSvgId: svgId
});

const features = [
  feature({
    resultKey: 'result-a', resultIndex: 0, resultName: 'a.svg',
    recordKey: 'record-a', biologicalFeatureId: 'feature-a', instanceHash: hashes.a,
    svgId: 'svg-a'
  }),
  feature({
    resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg',
    recordKey: 'record-b', biologicalFeatureId: 'feature-b', instanceHash: hashes.b,
    svgId: 'svg-b'
  }),
  feature({
    resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg',
    recordKey: 'record-b', biologicalFeatureId: 'feature-c', instanceHash: hashes.c,
    legendCaption: 'Accessory', svgId: 'svg-c'
  })
];

const targetKey = stableFeatureTargetKey(features[0]);
const selectedKey = stableFeatureTargetKey(features[2]);

const snapshot = (overrides = {}) => ({
  features: structuredClone(features),
  featureStrokeOverrides: {},
  legendStrokeOverrides: {},
  originalSvgStroke: { color: '#000000', width: 0.5 },
  selectedFeatureKeys: [targetKey, selectedKey],
  resultGenerationKey: 7,
  documentEpoch: 3,
  styleRevision: 4,
  styleFingerprint: 'sf1_fixture',
  exactScopeAvailable: true,
  ...overrides
});

const intent = (overrides = {}) => ({
  targetFeatureKey: targetKey,
  value: { kind: 'color', color: '#abcdef', strokeWidth: 2 },
  source: 'popup',
  originResultKey: 'result-a',
  originResultIndex: 0,
  resultGenerationKey: 7,
  documentEpoch: 3,
  styleRevision: 4,
  styleFingerprint: 'sf1_fixture',
  ...overrides
});

test('effective stroke resolves exact, caption, and SVG defaults per property', () => {
  const local = resolveFeatureStrokeViewModel({
    exactOverride: { strokeColor: '#112233' },
    captionOverride: { strokeColor: '#445566', strokeWidth: 1.5 },
    legendCaption: 'Core',
    svgDefaultStroke: { color: '#000000', width: 0.5 }
  });
  assert.deepEqual(local, {
    effectiveColor: '#112233',
    effectiveWidth: 1.5,
    explicitValue: '#112233',
    explicitWidth: null,
    origin: 'local',
    widthOrigin: 'legend-caption',
    originLabel: 'Feature override',
    canReset: true,
    allowNone: false,
    pickerEnabled: true
  });

  const inherited = resolveFeatureStrokeViewModel({
    captionOverride: { strokeColor: '#445566', strokeWidth: 1.5 },
    legendCaption: 'Core',
    svgDefaultStroke: { color: '#000000', width: 0.5 }
  });
  assert.equal(inherited.effectiveColor, '#445566');
  assert.equal(inherited.origin, 'legend-caption');
  assert.equal(inherited.explicitValue, null);
  assert.equal(inherited.canReset, false);
  assert.equal(inherited.pickerEnabled, true);
  assert.equal(inherited.allowNone, false);

  const fallback = resolveFeatureStrokeViewModel({
    svgDefaultStroke: { color: '#010203', width: 0.75 }
  });
  assert.equal(fallback.effectiveColor, '#010203');
  assert.equal(fallback.effectiveWidth, 0.75);
  assert.equal(fallback.origin, 'svg-default');
});

test('planner returns only legend-caption, selected-features, and single scopes', () => {
  const plan = planFeatureStrokeChange(snapshot(), intent());
  assert.equal(plan.status, 'needs-scope');
  assert.deepEqual(
    plan.candidates.map((candidate) => candidate.semanticScope),
    FEATURE_STROKE_SCOPE_ORDER
  );
  assert.equal(plan.candidates.some((candidate) => candidate.semanticScope === 'matching-rule'), false);

  const caption = plan.candidates[0];
  assert.equal(caption.resultExtent, 'all-affected-results');
  assert.equal(caption.targetCount, 2);
  assert.equal(caption.affectedResultCount, 2);
  assert.deepEqual(caption.targetsByResult.map((group) => ({
    resultKey: group.resultKey,
    count: group.featureKeys.length
  })), [{ resultKey: 'result-a', count: 1 }, { resultKey: 'result-b', count: 1 }]);

  const selected = plan.candidates[1];
  assert.equal(selected.resultExtent, 'all-affected-results');
  assert.deepEqual(selected.targetFeatureKeys, [targetKey, selectedKey].sort());
  assert.equal(plan.candidates[2].resultExtent, 'current-result');
  assert.match(caption.label, /2 features in 2 outputs/);
});

test('explicit group, selected, single, and inherit intents resolve immutably', () => {
  const pending = planFeatureStrokeChange(snapshot(), intent());
  const caption = pending.candidates.find((candidate) => candidate.semanticScope === 'legend-caption');
  const resolvedCaption = resolveFeatureStrokePlan(pending, caption.id);
  assert.equal(resolvedCaption.status, 'ready');
  assert.equal(resolvedCaption.semanticScope, 'legend-caption');
  assert.deepEqual(resolvedCaption.affectedResultKeys, ['result-a', 'result-b']);
  assert.equal(resolvedCaption.semanticAfter.durableStrokeIntent.kind, 'legend-caption');

  const selected = planFeatureStrokeChange(snapshot(), intent({
    source: 'selection-toolbar',
    targetFeatureKeys: [targetKey, selectedKey]
  }));
  assert.equal(selected.status, 'ready');
  assert.equal(resolveFeatureStrokePlan(selected).semanticScope, 'selected-features');

  const single = planFeatureStrokeChange(snapshot({ features: [features[0]] }), intent({
    semanticScope: 'single',
    value: { kind: 'inherit' }
  }));
  const resolvedSingle = resolveFeatureStrokePlan(single);
  assert.equal(resolvedSingle.semanticScope, 'single');
  assert.deepEqual(resolvedSingle.semanticAfter.value, { kind: 'inherit' });
  assert.equal(Object.isFrozen(resolvedSingle), true);
  assert.throws(() => {
    resolvedSingle.targetsByResult[0].featureKeys.push('mutated');
  }, TypeError);
  assert.equal(pending.status, 'needs-scope', 'resolution must not mutate the pending plan');
});

test('matching-rule is unavailable and stale revisions fail before candidate discovery', () => {
  const matching = planFeatureStrokeChange(snapshot(), intent({ semanticScope: 'matching-rule' }));
  assert.equal(matching.status, 'invalid');
  assert.equal(matching.diagnostics[0].code, 'unsupported-matching-rule');

  const stale = planFeatureStrokeChange(snapshot(), intent({ styleRevision: 5 }));
  assert.equal(stale.status, 'invalid');
  assert.equal(stale.diagnostics[0].code, 'stale-revision');
});

test('missing exact identity keeps group scope but fails closed for single and selection', () => {
  const legacyFeatures = structuredClone(features);
  delete legacyFeatures[0].instanceHash;
  delete legacyFeatures[2].instanceHash;
  const grouped = planFeatureStrokeChange(snapshot({
    features: legacyFeatures,
    exactScopeAvailable: false
  }), intent());
  assert.equal(grouped.status, 'needs-scope');
  assert.deepEqual(grouped.candidates.map((candidate) => candidate.semanticScope), ['legend-caption']);
  assert.deepEqual(grouped.diagnostics, [{
    code: 'exact-scope-unavailable',
    message: 'Generate to enable exact feature edits'
  }]);

  const isolated = planFeatureStrokeChange(snapshot({
    features: [legacyFeatures[2]],
    exactScopeAvailable: false,
    selectedFeatureKeys: []
  }), intent({
    targetFeatureKey: stableFeatureTargetKey(legacyFeatures[2]),
    originResultKey: 'result-b',
    originResultIndex: 1
  }));
  assert.equal(isolated.status, 'invalid');
  assert.equal(isolated.diagnostics[0].code, 'exact-scope-unavailable');
});

test('planning does not mutate deeply frozen inputs and is deterministic', () => {
  const frozenSnapshot = snapshot();
  Object.freeze(frozenSnapshot.features);
  frozenSnapshot.features.forEach(Object.freeze);
  Object.freeze(frozenSnapshot.featureStrokeOverrides);
  Object.freeze(frozenSnapshot.legendStrokeOverrides);
  Object.freeze(frozenSnapshot.originalSvgStroke);
  Object.freeze(frozenSnapshot.selectedFeatureKeys);
  Object.freeze(frozenSnapshot);
  const frozenIntent = Object.freeze(intent());
  assert.deepEqual(
    planFeatureStrokeChange(frozenSnapshot, frozenIntent),
    planFeatureStrokeChange(frozenSnapshot, frozenIntent)
  );
});
