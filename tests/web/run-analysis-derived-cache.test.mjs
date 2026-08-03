import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';

import {
  buildLosatDerivedPayloadCachePayload,
  hasRequiredCanonicalAnalysisResource
} from '../../gbdraw/web/js/app/run-analysis.js';

const baselineInput = {
  mode: 'collinear',
  maxHits: 5,
  bitscore: 50,
  evalue: '1e-5',
  identity: 70,
  alignmentLength: 0,
  collinearMinAnchors: 1,
  collinearMaxUnitGap: 0,
  collinearUnitMode: 'cds',
  collinearColorMode: 'orientation',
  collinearAnchorMode: 'rbh',
  collinearMaxDiagonalDrift: 0,
  collinearMaxConflictsInMergeGap: 1,
  collinearMaxParalogLinksPerOrthogroup: 2,
  collinearSearchScope: 'adjacent',
  orthogroupMembershipMode: 'anchor_core_v1',
  orthogroupMemberMaxHits: 5,
  recordPayloads: [
    {
      recordIndex: 1,
      proteinCacheKey: 'protein-b',
      runtimeBindingHash: 'runtime-b',
      displayBindingHash: 'display-b',
      viewTransform: { length: 200, reverse: true }
    },
    {
      recordIndex: 0,
      proteinCacheKey: 'protein-a',
      runtimeBindingHash: 'runtime-a',
      displayBindingHash: 'display-a',
      viewTransform: { length: 100, reverse: false }
    }
  ],
  pairPayloads: [
    {
      pairIndex: 0,
      queryIndex: 0,
      subjectIndex: 1,
      cacheKey: 'pair-a-b'
    }
  ]
};

const buildIdentity = (overrides = {}) => buildLosatDerivedPayloadCachePayload({
  ...baselineInput,
  ...overrides
});
const cacheKey = (identity) => createHash('sha256')
  .update(JSON.stringify(identity))
  .digest('hex');

const baselineIdentity = buildIdentity();
assert.equal(baselineIdentity.featureIdentity, 'stable-source-rendered-display-v1');
const unchangedIdentity = buildIdentity({
  recordPayloads: baselineInput.recordPayloads.map((record) => ({
    ...record,
    viewTransform: { ...record.viewTransform }
  })),
  pairPayloads: baselineInput.pairPayloads.map((pair) => ({ ...pair }))
});

assert.deepEqual(unchangedIdentity, baselineIdentity);
assert.equal(cacheKey(unchangedIdentity), cacheKey(baselineIdentity));

assert.equal(
  hasRequiredCanonicalAnalysisResource('collinear', {
    pairs: [],
    orthogroups: [],
    collinearGroups: []
  }),
  false,
  'derived cache payloads from before typed CollinearityResult propagation must miss'
);
assert.equal(
  hasRequiredCanonicalAnalysisResource('collinear', {
    collinearityResult: {
      schema: 2,
      kind: 'result',
      value: { type: 'CollinearityResult', fields: {} }
    }
  }),
  true
);
assert.equal(
  hasRequiredCanonicalAnalysisResource('orthogroup', {
    orthogroupResult: {
      schema: 1,
      kind: 'orthogroupResult',
      value: { type: 'OrthogroupResult', fields: {} }
    }
  }),
  true
);

const invalidationCases = [
  ['collinearMinAnchors', 'minAnchors', 2],
  ['collinearMaxUnitGap', 'maxGeneGap', 1],
  ['collinearMaxDiagonalDrift', 'maxDiagonalDrift', 1],
  ['collinearMaxConflictsInMergeGap', 'maxConflictsInMergeGap', 2],
  ['collinearMaxParalogLinksPerOrthogroup', 'maxParalogLinksPerOrthogroup', 3],
  ['collinearSearchScope', 'searchScope', 'all']
];

for (const [inputName, identityName, changedValue] of invalidationCases) {
  const changedIdentity = buildIdentity({ [inputName]: changedValue });
  assert.deepEqual(changedIdentity, {
    ...baselineIdentity,
    collinear: {
      ...baselineIdentity.collinear,
      [identityName]: String(changedValue)
    }
  });
  assert.notEqual(
    cacheKey(changedIdentity),
    cacheKey(baselineIdentity),
    `${inputName} must invalidate the derived cache key`
  );
}
