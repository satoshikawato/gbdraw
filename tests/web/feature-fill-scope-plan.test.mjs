import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

import {
  FEATURE_INSTANCE_HASH_QUALIFIER,
  normalizeFeatureFillValue,
  normalizePickerColor,
  normalizeRuleColor,
  resolveFeatureFillViewModel,
  resolveOrderedSpecificRule
} from '../../gbdraw/web/js/app/feature-editor/fill-view-model.js';
import {
  FEATURE_FILL_SCOPE_ORDER,
  planFeatureFillChange,
  resolveFeatureFillPlan,
  stableFeatureTargetKey,
  stableFeatureFillStringify
} from '../../gbdraw/web/js/app/feature-editor/fill-scope-plan.js';
import { parseSpecificRules } from '../../gbdraw/web/js/app/file-imports.js';

const hashes = {
  a: 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa',
  b: 'fi1_bbbbbbbbbbbbbbbbbbbbbbbbbb',
  c: 'fi1_cccccccccccccccccccccccccc',
  d: 'fi1_dddddddddddddddddddddddddd',
  e: 'fi1_eeeeeeeeeeeeeeeeeeeeeeeeee'
};

const precedenceOracle = JSON.parse(await readFile(
  fileURLToPath(new URL('../fixtures/specific_color_precedence_oracle.json', import.meta.url)),
  'utf8'
));

const feature = ({
  resultKey,
  resultName,
  resultIndex,
  recordKey,
  biologicalFeatureId,
  instanceHash,
  type = 'CDS',
  product,
  renderedLabel,
  legendCaption,
  orthogroupId = '',
  rendered = true
}) => ({
  resultKey,
  resultName,
  resultIndex,
  recordKey,
  biologicalFeatureId,
  instanceHash,
  type,
  product,
  qualifiers: product ? { product: [product] } : {},
  renderedLabel,
  legendCaption,
  orthogroupId,
  ...(rendered ? { renderedFeatureSvgId: `svg-${biologicalFeatureId}` } : { rendered: false })
});

const features = [{
  resultKey: 'result-a', resultName: 'a.svg', resultIndex: 0,
  recordKey: 'record-a', biologicalFeatureId: 'feature-a', instanceHash: hashes.a,
  product: 'Kinase', renderedLabel: 'Shared rendered label', legendCaption: 'Core',
  orthogroupId: 'og-shared'
}, {
  resultKey: 'result-b', resultName: 'b.svg', resultIndex: 1,
  recordKey: 'record-b', biologicalFeatureId: 'feature-b', instanceHash: hashes.b,
  product: 'Kinase', renderedLabel: 'Other label', legendCaption: 'Core'
}, {
  resultKey: 'result-c', resultName: 'c.svg', resultIndex: 2,
  recordKey: 'record-c', biologicalFeatureId: 'feature-c', instanceHash: hashes.c,
  product: 'Kinase', renderedLabel: '', legendCaption: 'Core', rendered: false
}, {
  resultKey: 'result-b', resultName: 'b.svg', resultIndex: 1,
  recordKey: 'record-b', biologicalFeatureId: 'feature-d', instanceHash: hashes.d,
  type: 'tRNA', product: 'Transfer RNA', renderedLabel: 'Shared rendered label', legendCaption: 'tRNA',
  orthogroupId: 'og-shared'
}, {
  resultKey: 'result-a', resultName: 'a.svg', resultIndex: 0,
  recordKey: 'record-a', biologicalFeatureId: 'feature-e', instanceHash: hashes.e,
  type: 'misc_feature', product: 'Marker', renderedLabel: 'Marker', legendCaption: 'Marker'
}].map(feature);

const targetAlias = 'record-a\u0000feature-a';
const rules = [{
  feat: 'CDS', qual: 'product', val: '^Kinase$', color: '#123456', cap: 'Core'
}];

const baseSnapshot = () => ({
  features: structuredClone(features),
  specificRules: structuredClone(rules),
  paletteColors: { CDS: '#cccccc', tRNA: '#fedcba', misc_feature: '#334455' },
  selectedFeatureKeys: [targetAlias, 'record-a\u0000feature-e'],
  resultGenerationKey: 7,
  documentEpoch: 3,
  styleFingerprint: 'style-1',
  exactScopeAvailable: true
});

const intent = (overrides = {}) => ({
  targetFeatureKey: targetAlias,
  value: { kind: 'color', color: '#abcdef' },
  requestedCaption: 'Updated Core',
  source: 'popup',
  originResultKey: 'result-a',
  originResultIndex: 0,
  resultGenerationKey: 7,
  documentEpoch: 3,
  styleFingerprint: 'style-1',
  ...overrides
});

const deepFreeze = (value) => {
  if (!value || typeof value !== 'object' || Object.isFrozen(value)) return value;
  Object.values(value).forEach(deepFreeze);
  return Object.freeze(value);
};

test('effective fill view model keeps inherited colors editable and reports origin separately', () => {
  assert.deepEqual(resolveFeatureFillViewModel({
    feature: { type: 'CDS' },
    paletteColors: { CDS: '#54BCF8' }
  }), {
    effectiveColor: '#54bcf8',
    explicitValue: null,
    origin: 'palette',
    originLabel: 'Inherited from palette',
    canReset: false,
    allowNone: true
  });

  assert.deepEqual(resolveFeatureFillViewModel({
    feature: features[0],
    specificRules: rules,
    paletteColors: { CDS: '#cccccc' }
  }), {
    effectiveColor: '#123456',
    explicitValue: null,
    origin: 'specific-rule',
    originLabel: 'Inherited from rule: Core',
    canReset: false,
    allowNone: true
  });

  assert.deepEqual(resolveFeatureFillViewModel({
    feature: { type: 'CDS' },
    explicitValue: '#ff00ff',
    explicitOrigin: 'derived-replay',
    paletteColors: { CDS: '#54bcf8' }
  }), {
    effectiveColor: '#54bcf8',
    explicitValue: null,
    origin: 'palette',
    originLabel: 'Inherited from palette',
    canReset: false,
    allowNone: true
  });

  assert.deepEqual(resolveFeatureFillViewModel({
    feature: { type: 'CDS' },
    localRule: { feat: 'CDS', qual: 'product', val: '.*', color: '#ff00ff' },
    paletteColors: { CDS: '#54bcf8' }
  }), {
    effectiveColor: '#54bcf8',
    explicitValue: null,
    origin: 'palette',
    originLabel: 'Inherited from palette',
    canReset: false,
    allowNone: true
  });

  const exact = {
    feat: 'CDS', qual: FEATURE_INSTANCE_HASH_QUALIFIER,
    val: hashes.a, color: 'none', cap: 'Hidden kinase'
  };
  assert.deepEqual(resolveFeatureFillViewModel({
    feature: features[0],
    specificRules: [exact, ...rules],
    paletteColors: { CDS: '#cccccc' },
    allowNone: false
  }), {
    effectiveColor: 'none',
    explicitValue: 'none',
    origin: 'local',
    originLabel: 'Feature override',
    canReset: true,
    allowNone: false
  });
});

test('fill values preserve inherit, none, legacy RuleColor forms, and picker domain', () => {
  assert.deepEqual(normalizeFeatureFillValue(null), { kind: 'inherit' });
  assert.deepEqual(normalizeFeatureFillValue({ kind: 'inherit' }), { kind: 'inherit' });
  assert.deepEqual(normalizeFeatureFillValue('none'), { kind: 'none' });
  for (const color of ['#abc', '#abcd', '#abcdef', '#abcdef12']) {
    assert.equal(normalizeRuleColor(color.toUpperCase()), color);
    assert.deepEqual(normalizeFeatureFillValue({ kind: 'color', color }), { kind: 'color', color });
  }
  assert.equal(normalizeRuleColor('#abcde'), null);
  assert.equal(normalizePickerColor('#abcdef'), '#abcdef');
  assert.equal(normalizePickerColor('#abc'), null);
  assert.equal(normalizePickerColor('none'), null);
});

test('specific-rule precedence nests exact type before wildcard and reserves the exact literal', () => {
  const target = {
    ...features[0],
    hash: 'public-hash',
    qualifiers: { ...features[0].qualifiers, instance_hash: ['biological-value'] }
  };
  const exactTypeHash = { feat: 'CDS', qual: 'hash', val: '^public-hash$', color: '#222222' };
  const wildcardInstance = {
    feat: '*', qual: FEATURE_INSTANCE_HASH_QUALIFIER, val: hashes.a, color: '#111111'
  };
  assert.equal(
    resolveOrderedSpecificRule(target, [wildcardInstance, exactTypeHash]).rule,
    exactTypeHash
  );

  const exactTypeInstance = {
    feat: 'CDS', qual: FEATURE_INSTANCE_HASH_QUALIFIER, val: hashes.a, color: '#333333'
  };
  assert.equal(
    resolveOrderedSpecificRule(target, [exactTypeHash, exactTypeInstance]).rule,
    exactTypeInstance
  );

  const biologicalQualifier = {
    feat: 'CDS', qual: 'instance_hash', val: '^biological-', color: '#444444'
  };
  assert.equal(resolveOrderedSpecificRule(target, [biologicalQualifier]).rule, biologicalQualifier);
  assert.equal(
    resolveOrderedSpecificRule(target, [{
      ...exactTypeInstance,
      qual: '__GBDRAW_INSTANCE_HASH__'
    }]),
    null
  );
});

test('schema-5 reserved semantic spelling remains a biological regex qualifier', () => {
  const parsed = parseSpecificRules(
    'CDS\t__gbdraw_semantic_scope__\tbio.*\t#abcdef\tBiological value\n',
    { schema: 5 }
  ).rules;
  assert.equal(parsed[0].match, 'regex');
  assert.equal(resolveOrderedSpecificRule({
    type: 'CDS',
    qualifiers: { __gbdraw_semantic_scope__: ['bio-value'] }
  }, parsed)?.rule, parsed[0]);
});

test('Web specific-color resolution matches the shared Python precedence oracle', () => {
  const feature = {
    type: precedenceOracle.feature.type,
    record_id: precedenceOracle.feature.recordId,
    recordKey: precedenceOracle.feature.recordKey,
    start: precedenceOracle.feature.start,
    end: precedenceOracle.feature.end,
    strand: precedenceOracle.feature.strand,
    hash: precedenceOracle.feature.publicHash,
    stableFeatureId: precedenceOracle.feature.publicHash,
    instanceHash: precedenceOracle.feature.instanceHash,
    qualifiers: precedenceOracle.feature.qualifiers
  };
  precedenceOracle.cases.forEach((entry) => {
    const rules = entry.rules.map(([feat, qual, val, color, cap, match]) => ({
      feat, qual, val, color, cap, match
    }));
    assert.equal(
      resolveOrderedSpecificRule(feature, rules)?.rule?.cap,
      entry.expectedCaption,
      entry.id
    );
  });
});

test('planner returns complete deterministic scope order and batch-owned targets', () => {
  const snapshot = deepFreeze(baseSnapshot());
  const request = deepFreeze(intent());
  const beforeSnapshot = structuredClone(snapshot);
  const beforeIntent = structuredClone(request);
  const first = planFeatureFillChange(snapshot, request);
  const second = planFeatureFillChange(snapshot, request);

  assert.equal(first.status, 'needs-scope');
  assert.deepEqual(first, second);
  assert.deepEqual(
    first.candidates.map((candidate) => candidate.semanticScope),
    FEATURE_FILL_SCOPE_ORDER
  );
  assert.deepEqual(snapshot, beforeSnapshot);
  assert.deepEqual(request, beforeIntent);
  assert.doesNotThrow(() => JSON.parse(JSON.stringify(first)));

  const matching = first.candidates[0];
  assert.equal(matching.resultExtent, 'all-affected-results');
  assert.equal(matching.targetCount, 3);
  assert.equal(matching.affectedResultCount, 3);
  assert.deepEqual(
    matching.targetsByResult.map(({ resultKey, featureKeys }) => [resultKey, featureKeys.length]),
    [['result-a', 1], ['result-b', 1], ['result-c', 1]]
  );
  assert.equal(matching.durableRuleIntent.kind, 'update-specific-rule');

  const featureType = first.candidates.find(
    (candidate) => candidate.semanticScope === 'feature-type'
  );
  assert.equal(featureType.targetCount, 3);
  assert.equal(featureType.affectedResultCount, 3);
  assert.deepEqual(featureType.durableRuleIntent, {
    kind: 'feature-type', featureType: 'CDS'
  });

  const similarityGroup = first.candidates.find(
    (candidate) => candidate.semanticScope === 'similarity-group'
  );
  assert.equal(similarityGroup.targetCount, 2);
  assert.equal(similarityGroup.affectedResultCount, 2);
  assert.deepEqual(similarityGroup.durableRuleIntent, {
    kind: 'similarity-group', groupId: 'og-shared'
  });

  const selected = first.candidates.find(
    (candidate) => candidate.semanticScope === 'selected-features'
  );
  assert.equal(selected.resultExtent, 'current-result');
  assert.equal(selected.targetCount, 2);
  assert.equal(selected.durableRuleIntent.selectors.every(
    (selector) => selector.qual === FEATURE_INSTANCE_HASH_QUALIFIER
  ), true);
  assert.equal(first.candidates.at(-1).targetCount, 1);
});

test('Similarity candidates preserve every target multi-membership as distinct durable intent', () => {
  const target = {
    ...features[0],
    orthogroupId: '',
    orthogroupIds: ['og-alpha', 'og-beta']
  };
  const alphaPeer = {
    ...features[1],
    orthogroupId: 'og-alpha'
  };
  const betaPeer = {
    ...features[3],
    orthogroupId: 'og-beta'
  };
  const snapshot = {
    ...baseSnapshot(),
    features: [target, alphaPeer, betaPeer]
  };
  const plan = planFeatureFillChange(snapshot, intent({
    targetFeatureKey: stableFeatureTargetKey(target)
  }));
  const similarity = plan.candidates.filter(
    (candidate) => candidate.semanticScope === 'similarity-group'
  );

  assert.deepEqual(
    similarity.map((candidate) => candidate.durableRuleIntent.groupId),
    ['og-alpha', 'og-beta']
  );
  assert.equal(new Set(similarity.map((candidate) => candidate.id)).size, 2);
  assert.deepEqual(similarity.map((candidate) => candidate.targetCount), [2, 2]);
});

test('planner derives rendered-label and Similarity-group scopes from production catalog owners', () => {
  const item = {
    resultKey: 'catalog-result',
    resultName: 'catalog.svg',
    biologicalFeatures: [{
      recordKey: 'record-a', biologicalFeatureId: 'feature-a', instanceHash: hashes.a,
      type: 'CDS', qualifiers: { product: ['Kinase'] }
    }, {
      recordKey: 'record-b', biologicalFeatureId: 'feature-b', instanceHash: hashes.b,
      type: 'CDS', qualifiers: { product: ['Kinase'] }
    }, {
      recordKey: 'record-c', biologicalFeatureId: 'feature-c', instanceHash: hashes.c,
      type: 'CDS', qualifiers: { product: ['Other'] }
    }],
    features: [{
      recordKey: 'record-a', biologicalFeatureId: 'feature-a', svgId: 'svg-a'
    }, {
      recordKey: 'record-b', biologicalFeatureId: 'feature-b', svgId: 'svg-b'
    }, {
      recordKey: 'record-c', biologicalFeatureId: 'feature-c', svgId: 'svg-c'
    }],
    orthogroups: [{
      id: 'og-catalog', presentationScope: 'global_collinear', groupKind: 'orthogroup',
      members: [
        { recordKey: 'record-a', biologicalFeatureId: 'feature-a' },
        { recordKey: 'record-b', biologicalFeatureId: 'feature-b' }
      ]
    }]
  };
  const withResult = (entry, index) => ({
    ...entry,
    ...item.features[index],
    resultKey: item.resultKey,
    resultName: item.resultName,
    resultIndex: 0,
    renderedFeatureSvgId: item.features[index].svgId
  });
  const snapshot = {
    catalog: { schema: 3, items: [item] },
    features: item.biologicalFeatures.map((entry, index) => withResult(entry, index)),
    extractedFeatures: item.biologicalFeatures.map((entry, index) => withResult(entry, index)),
    editableLabels: [
      { featureId: 'svg-a', text: 'Shared rendered label' },
      { featureId: 'svg-b', text: 'Shared rendered label' }
    ],
    specificRules: [],
    paletteColors: { CDS: '#cccccc' },
    exactScopeAvailable: true
  };
  const target = stableFeatureTargetKey(snapshot.features[0]);
  const plan = planFeatureFillChange(snapshot, {
    targetFeatureKey: target,
    value: { kind: 'color', color: '#abcdef' },
    source: 'popup',
    originResultKey: item.resultKey,
    originResultIndex: 0
  });

  assert.ok(['needs-scope', 'conflict'].includes(plan.status));
  const candidates = Object.fromEntries(
    plan.candidates.map((candidate) => [candidate.semanticScope, candidate])
  );
  assert.equal(candidates['feature-type'].targetCount, 3);
  assert.equal(candidates['rendered-label'].targetCount, 2);
  assert.equal(candidates['source-annotation-label'].targetCount, 2);
  assert.equal(candidates['similarity-group'].targetCount, 2);
  const pendingConflict = resolveFeatureFillPlan(plan, candidates['similarity-group'].id);
  assert.equal(pendingConflict.status, 'conflict');
  assert.equal(pendingConflict.selectedCandidateId, candidates['similarity-group'].id);
  assert.equal(resolveFeatureFillPlan(plan, candidates['similarity-group'].id, {
    conflictChoiceId: 'separate-caption'
  }).status, 'ready');
});

test('source-annotation scope groups equal labels across different source qualifiers', () => {
  const snapshot = baseSnapshot();
  snapshot.features = [
    feature({
      resultKey: 'result-a', resultName: 'a.svg', resultIndex: 0,
      recordKey: 'record-a', biologicalFeatureId: 'feature-a', instanceHash: hashes.a,
      product: 'Shared annotation', renderedLabel: 'Product label', legendCaption: 'Core'
    }),
    {
      ...feature({
        resultKey: 'result-b', resultName: 'b.svg', resultIndex: 1,
        recordKey: 'record-b', biologicalFeatureId: 'feature-b', instanceHash: hashes.b,
        renderedLabel: 'Gene label', legendCaption: 'Core'
      }),
      gene: 'Shared annotation',
      qualifiers: { gene: ['Shared annotation'] }
    }
  ];
  snapshot.specificRules = [];
  snapshot.selectedFeatureKeys = [];

  const plan = planFeatureFillChange(snapshot, intent({
    requestedCaption: 'Core',
    targetFeatureKey: 'record-a\u0000feature-a'
  }));
  const candidate = plan.candidates.find(
    (entry) => entry.semanticScope === 'source-annotation-label'
  );
  assert.equal(candidate.targetCount, 2);
  assert.equal(candidate.affectedResultCount, 2);
});

test('group scopes keep semantic extent when only one Result happens to be affected', () => {
  const snapshot = baseSnapshot();
  snapshot.features = snapshot.features.filter((entry) => entry.resultKey === 'result-a');
  snapshot.features[1].legendCaption = 'Core';
  snapshot.selectedFeatureKeys = [];
  const plan = planFeatureFillChange(snapshot, intent());
  const caption = plan.candidates.find((candidate) => candidate.semanticScope === 'legend-caption');

  assert.equal(caption.resultExtent, 'all-affected-results');
  assert.equal(caption.affectedResultCount, 1);
  assert.equal(caption.targetCount, 2);
});

test('planner auto-resolves single only when no group exists and exact identity is valid', () => {
  const snapshot = baseSnapshot();
  snapshot.features = [snapshot.features[4]];
  snapshot.selectedFeatureKeys = [];
  const request = intent({
    targetFeatureKey: 'record-a\u0000feature-e',
    originResultKey: 'result-a',
    requestedCaption: 'Marker changed'
  });
  const plan = planFeatureFillChange(snapshot, request);

  assert.equal(plan.status, 'ready');
  assert.deepEqual(plan.candidates.map((candidate) => candidate.semanticScope), ['single']);
  assert.equal(plan.selectedCandidateId, plan.candidates[0].id);

  const resolved = resolveFeatureFillPlan(plan);
  assert.equal(resolved.status, 'ready');
  assert.equal(resolved.semanticScope, 'single');
  assert.deepEqual(resolved.affectedResultKeys, ['result-a']);
  assert.equal(resolved.semanticAfter.durableRuleIntent.selectors[0].val, hashes.e);
});

test('legacy catalogues retain group planning but fail closed for exact scopes', () => {
  const snapshot = baseSnapshot();
  snapshot.features.forEach((entry) => delete entry.instanceHash);
  snapshot.exactScopeAvailable = false;
  const plan = planFeatureFillChange(snapshot, intent());

  assert.equal(plan.status, 'needs-scope');
  assert.equal(plan.candidates.some((candidate) => candidate.semanticScope === 'single'), false);
  assert.equal(
    plan.candidates.some((candidate) => candidate.semanticScope === 'selected-features'),
    false
  );
  assert.deepEqual(plan.diagnostics, [{
    code: 'exact-scope-unavailable',
    message: 'Generate to enable exact feature edits'
  }]);

  snapshot.features = [snapshot.features[4]];
  const unavailable = planFeatureFillChange(snapshot, intent({
    targetFeatureKey: 'record-a\u0000feature-e', requestedCaption: 'Marker changed'
  }));
  assert.equal(unavailable.status, 'invalid');
});

test('stale identity and broad current-Result snapshots are rejected', () => {
  assert.equal(
    planFeatureFillChange(baseSnapshot(), intent({ resultGenerationKey: 8 })).diagnostics[0].code,
    'stale-intent'
  );
  assert.equal(
    planFeatureFillChange(baseSnapshot(), intent({ originResultKey: 'result-b' })).diagnostics[0].code,
    'stale-result'
  );
  assert.equal(
    planFeatureFillChange(baseSnapshot(), intent({
      semanticScope: 'legend-caption', resultExtent: 'current-result'
    })).diagnostics[0].code,
    'unsupported-current-result-snapshot'
  );
});

test('caption color conflicts are explicit and resolve without changing scope', () => {
  const snapshot = baseSnapshot();
  snapshot.features.push(feature({
    resultKey: 'result-b', resultName: 'b.svg', resultIndex: 1,
    recordKey: 'record-b', biologicalFeatureId: 'existing', instanceHash: hashes.d,
    type: 'gene', product: 'Existing', renderedLabel: 'Existing', legendCaption: 'Collision'
  }));
  snapshot.features.at(-1).explicitFillValue = '#112233';
  const plan = planFeatureFillChange(snapshot, intent({ requestedCaption: 'Collision' }));
  assert.equal(plan.status, 'conflict');
  assert.deepEqual(plan.conflict.choices.map((choice) => choice.id), [
    'reuse-existing', 'separate-caption'
  ]);

  const selected = plan.candidates.find((candidate) => candidate.semanticScope === 'single');
  const reused = resolveFeatureFillPlan(plan, selected.id, { conflictChoiceId: 'reuse-existing' });
  assert.equal(reused.status, 'ready');
  assert.deepEqual(reused.intent.value, { kind: 'color', color: '#112233' });
  assert.equal(reused.semanticScope, 'single');

  const separate = resolveFeatureFillPlan(plan, selected.id, { conflictChoiceId: 'separate-caption' });
  assert.equal(separate.intent.requestedCaption, 'Collision (2)');
  assert.equal(plan.status, 'conflict', 'resolution must not mutate the pending plan');
});

test('blank UI captions inherit the current caption and require conflict choice only for exact scopes', () => {
  const plan = planFeatureFillChange(baseSnapshot(), intent({ requestedCaption: '' }));
  assert.equal(plan.status, 'conflict');
  assert.equal(plan.intent.requestedCaption, 'Core');

  const group = plan.candidates.find((candidate) => candidate.semanticScope === 'legend-caption');
  const grouped = resolveFeatureFillPlan(plan, group.id);
  assert.equal(grouped.status, 'ready');
  assert.equal(grouped.semanticScope, 'legend-caption');
  assert.equal(grouped.semanticAfter.requestedCaption, 'Core');

  const single = plan.candidates.find((candidate) => candidate.semanticScope === 'single');
  assert.equal(resolveFeatureFillPlan(plan, single.id).status, 'conflict');
  const separated = resolveFeatureFillPlan(plan, single.id, {
    conflictChoiceId: 'separate-caption'
  });
  assert.equal(separated.status, 'ready');
  assert.equal(separated.semanticAfter.requestedCaption, 'Core (2)');
});

test('stable target keys do not persist navigation index as selector identity', () => {
  const original = feature({
    resultKey: 'durable-result', resultName: 'first.svg', resultIndex: 0,
    recordKey: 'record', biologicalFeatureId: 'feature', instanceHash: hashes.a,
    product: 'A', renderedLabel: 'A', legendCaption: 'A'
  });
  const moved = { ...original, resultIndex: 9, resultName: 'renamed.svg' };
  assert.equal(stableFeatureTargetKey(original), stableFeatureTargetKey(moved));
  assert.equal(
    stableFeatureTargetKey(original),
    stableFeatureFillStringify(['durable-result', 'record', 'feature'])
  );
});
