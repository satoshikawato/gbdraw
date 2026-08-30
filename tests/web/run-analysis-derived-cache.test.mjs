import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

import {
  buildLosatDerivedPayloadCachePayload,
  hasRequiredCanonicalAnalysisResource,
  resolveProteinBlastpCandidateLimit
} from '../../gbdraw/web/js/app/run-analysis.js';

assert.equal(resolveProteinBlastpCandidateLimit(5), 5);
assert.equal(resolveProteinBlastpCandidateLimit('7'), 7);
assert.equal(resolveProteinBlastpCandidateLimit(null), null);
assert.equal(resolveProteinBlastpCandidateLimit(undefined), null);
assert.throws(() => resolveProteinBlastpCandidateLimit(0), /Candidate limit/);
assert.throws(() => resolveProteinBlastpCandidateLimit('invalid'), /Candidate limit/);

const runAnalysisPath = resolve('gbdraw/web/js/app/run-analysis.js');
const runAnalysisUrl = pathToFileURL(runAnalysisPath);
const identityProbeDir = await mkdtemp(join(tmpdir(), 'gbdraw-raw-identity-'));
const identityProbePath = join(identityProbeDir, 'run-analysis-identity-probe.mjs');
const identityProbeSource = (await readFile(runAnalysisPath, 'utf8'))
  .replace(
    'const buildLosatCachePayload = ({',
    'export const buildLosatCachePayload = ({'
  )
  .replace(/from '(\.\.?\/[^']+)'/g, (_match, specifier) => (
    `from '${new URL(specifier, runAnalysisUrl).href}'`
  ));
await writeFile(identityProbePath, identityProbeSource, 'utf8');
const { buildLosatCachePayload } = await import(pathToFileURL(identityProbePath));

const rawIdentityInput = {
  identityKind: 'protein',
  program: 'blastp',
  outfmt: '6',
  args: ['--max-target-seqs', '9'],
  queryProteinSetHash: 'query-proteins',
  subjectProteinSetHash: 'subject-proteins',
  queryRuntimeBindingHash: 'query-runtime',
  subjectRuntimeBindingHash: 'subject-runtime',
  queryRecordInstanceKey: 'query-record',
  subjectRecordInstanceKey: 'subject-record'
};
const rawIdentity = buildLosatCachePayload(rawIdentityInput);
for (const appearanceOnly of [
  { pairwiseMatchStyle: 'curve' },
  { comparisonHeight: 90 },
  { comparisonDisclosureOpen: true },
  { scaleFontSize: 27 },
  { rulerLabelFontSize: 11 }
]) {
  assert.deepEqual(
    buildLosatCachePayload({ ...rawIdentityInput, ...appearanceOnly }),
    rawIdentity,
    'render appearance and disclosure state must not change raw scientific identity'
  );
}

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
  collinearMergeOrientation: 'either',
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
  ['orthogroupMemberMaxHits', 'memberMaxHits', 7, 'orthogroup'],
  ['collinearMinAnchors', 'minAnchors', 2],
  ['collinearMaxUnitGap', 'maxGeneGap', 1],
  ['collinearUnitMode', 'unitMode', 'locus'],
  ['collinearColorMode', 'colorMode', 'orientation_identity'],
  ['collinearAnchorMode', 'anchorMode', 'all'],
  ['collinearMergeOrientation', 'mergeOrientation', 'strand'],
  ['collinearMaxDiagonalDrift', 'maxDiagonalDrift', 1],
  ['collinearMaxConflictsInMergeGap', 'maxConflictsInMergeGap', 2],
  ['collinearMaxParalogLinksPerOrthogroup', 'maxParalogLinksPerOrthogroup', 3],
  ['collinearSearchScope', 'searchScope', 'all']
];

for (const [inputName, identityName, changedValue, section = 'collinear'] of invalidationCases) {
  const changedIdentity = buildIdentity({ [inputName]: changedValue });
  assert.deepEqual(changedIdentity, {
    ...baselineIdentity,
    [section]: {
      ...baselineIdentity[section],
      [identityName]: section === 'orthogroup' ? Number(changedValue) : String(changedValue)
    }
  });
  assert.notEqual(
    cacheKey(changedIdentity),
    cacheKey(baselineIdentity),
    `${inputName} must invalidate the derived cache key`
  );
  assert.deepEqual(changedIdentity.pairs, baselineIdentity.pairs);
}

assert.deepEqual(
  buildIdentity({ maxHits: 99, pairwiseMatchStyle: 'curve' }),
  baselineIdentity,
  'Collinear derived identity must ignore Pairwise-only and match-style settings'
);
for (const appearanceOnly of [
  { pairwiseMatchStyle: 'curve' },
  { comparisonHeight: 90 },
  { comparisonDisclosureOpen: true },
  { scaleFontSize: 27 },
  { rulerLabelFontSize: 11 }
]) {
  assert.deepEqual(
    buildIdentity(appearanceOnly),
    baselineIdentity,
    'render appearance and disclosure state must not change derived scientific identity'
  );
}

const orthogroupIdentity = buildIdentity({ mode: 'orthogroup' });
assert.equal(Object.hasOwn(orthogroupIdentity, 'pairwise'), false);
assert.equal(Object.hasOwn(orthogroupIdentity, 'collinear'), false);
assert.deepEqual(
  buildIdentity({
    mode: 'orthogroup',
    maxHits: 99,
    collinearUnitMode: 'locus',
    collinearAnchorMode: 'all',
    collinearMergeOrientation: 'strand',
    collinearSearchScope: 'all'
  }),
  orthogroupIdentity,
  'Orthogroup derived identity must ignore Pairwise-only and Collinear-only settings'
);

const pairwiseIdentity = buildIdentity({ mode: 'pairwise' });
assert.equal(Object.hasOwn(pairwiseIdentity, 'orthogroup'), false);
assert.equal(Object.hasOwn(pairwiseIdentity, 'collinear'), false);
assert.notEqual(
  cacheKey(buildIdentity({ mode: 'pairwise', maxHits: 6 })),
  cacheKey(pairwiseIdentity),
  'Pairwise max hits must participate in Pairwise derived identity'
);
