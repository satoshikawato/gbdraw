import assert from 'node:assert/strict';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

import * as runAnalysis from '../../gbdraw/web/js/app/run-analysis.js';
import { bindCanonicalTypedResource } from '../../gbdraw/web/js/services/session-request.js';

const { hasRequiredCanonicalAnalysisResource, resolveProteinBlastpCandidateLimit } = runAnalysis;

assert.equal(resolveProteinBlastpCandidateLimit(5), 5);
assert.equal(resolveProteinBlastpCandidateLimit('7'), 7);
assert.equal(resolveProteinBlastpCandidateLimit(null), null);
assert.equal(resolveProteinBlastpCandidateLimit(undefined), null);
assert.throws(() => resolveProteinBlastpCandidateLimit(0), /Max target seqs/);
assert.throws(() => resolveProteinBlastpCandidateLimit('invalid'), /Max target seqs/);

const runAnalysisPath = resolve('gbdraw/web/js/app/run-analysis.js');
const runAnalysisUrl = pathToFileURL(runAnalysisPath);
const identityProbeDir = await mkdtemp(join(tmpdir(), 'gbdraw-raw-identity-'));
const identityProbePath = join(identityProbeDir, 'run-analysis-identity-probe.mjs');
const identityProbeSource = (await readFile(runAnalysisPath, 'utf8'))
  .replace('const stripRuntimeCacheStats =', 'export const stripRuntimeCacheStats =')
  .replace(/from '(\.\.?\/[^']+)'/g, (_match, specifier) => (
    `from '${new URL(specifier, runAnalysisUrl).href}'`
  ));
await writeFile(identityProbePath, identityProbeSource, 'utf8');
const { buildLosatCachePayload, stripRuntimeCacheStats } = await import(pathToFileURL(identityProbePath));

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

// OV-1063: the derived identity and key have one owner, gbdraw.analysis.protein_artifacts
// (tests/test_derived_losatp_identity.py). Generate looks a payload up by the derived key of the
// `buildProteinLosatCacheKeys` reply, and the Web hashes no derived identity of its own.
assert.equal(Object.hasOwn(runAnalysis, 'buildLosatDerivedPayloadCachePayload'), false);
{
  const source = await readFile(runAnalysisPath, 'utf8');
  assert.match(source, /\(\{ keys: proteinCacheKeys, derivedKey: derivedCacheKey \} = requireProteinLosatCacheKeys\(/);
  assert.match(source, /getLosatDerivedCacheEntry\(\s*derivedCacheMap,\s*derivedCacheKey,/);
  assert.doesNotMatch(source, /derived-option-conformance/);
}

// Review F3: a converter whose identity key differs from the looked-up key is
// a classified Result failure that names its reason. The check runs as written
// in run-analysis.js, on the converter's reply.
{
  const { diagnosticError, normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');
  const source = await readFile(runAnalysisPath, 'utf8');
  const start = source.indexOf('              convertedPayload = /** @type {Record<string, any>} */ (response.result);');
  const end = source.indexOf('              if (["orthogroup", "collinear"].includes(blastpMode)) {', start);
  assert.ok(start > 0 && end > start);
  const admit = new Function(
    'response', 'useDerivedProteinPayloadCache', 'derivedCacheKey', 'diagnosticError',
    `let convertedPayload;\n${source.slice(start, end)}\nreturn convertedPayload;`
  );
  assert.deepEqual(admit({ result: { pairs: [], identityKey: 'a'.repeat(64) } }, true, 'a'.repeat(64), diagnosticError), { pairs: [] });
  let thrown = null;
  try {
    admit({ result: { pairs: [], identityKey: 'b'.repeat(64) } }, true, 'a'.repeat(64), diagnosticError);
  } catch (error) {
    thrown = error;
  }
  const failure = normalizeUserFacingError(thrown);
  assert.equal(failure.code, 'RESULT_INVALID');
  assert.equal(failure.stage, 'result-admission');
  assert.equal(failure.operation, 'convertLosatpPairsToGenomicPayload');
  assert.equal(failure.context.reason, 'DERIVED_KEY_MISMATCH');
  assert.match(failure.summary, /does not match the settings it was looked up with/);
}

assert.equal(
  hasRequiredCanonicalAnalysisResource('collinear', {
    pairs: [],
    orthogroups: [],
    collinearGroups: []
  }),
  false,
  'derived cache payloads from before typed CollinearityResult propagation must miss'
);
// The typed result as the analysis helper returns it: parsed from its Python bytes and bound to them.
const boundTypedResult = (typed, kind) => {
  const text = JSON.stringify(typed);
  bindCanonicalTypedResource(typed, {
    kind, type: 'application/json', size: Buffer.byteLength(text), lastModified: 0,
    encoding: 'base64', data: Buffer.from(text).toString('base64')
  });
  return typed;
};
assert.equal(
  hasRequiredCanonicalAnalysisResource('collinear', {
    pairs: [],
    collinearityResult: boundTypedResult({
      schema: 3,
      kind: 'result',
      value: { type: 'CollinearityResult', fields: {} }
    }, 'collinearity-result')
  }),
  true
);
assert.equal(
  hasRequiredCanonicalAnalysisResource('orthogroup', {
    pairs: [],
    orthogroupResult: boundTypedResult({
      schema: 3,
      kind: 'orthogroupResult',
      value: { type: 'OrthogroupGraphResult', fields: {} }
    }, 'orthogroup-result')
  }),
  true
);

// Old persisted resources remain readable, but do not masquerade as new cache hits.
for (const schema of [1, 2]) {
  assert.equal(hasRequiredCanonicalAnalysisResource('orthogroup', {
    orthogroupResult: {schema, kind: 'orthogroupResult', value: {type: 'OrthogroupResult', fields: {}}}
  }), false);
}
// OV-1040 / D-66: a derived payload is a hit only when its typed result is bound to the Python bytes and every
// pair holds its canonical table. A result read from a Session is not bound, and a pair with only `rows` predates
// the canonical table: both are misses, so the analysis runs again instead of encoding in JS.
{
  const typedText = '{"kind":"orthogroupResult","schema":3,"value":{"fields":{"score":1.0},"type":"OrthogroupGraphResult"}}';
  const orthogroupResult = JSON.parse(typedText);
  const payload = { pairs: [{ pair_index: 0, tsv: '', canonicalTsv: 'query\n', hit_count: 0 }], orthogroupResult };
  assert.equal(hasRequiredCanonicalAnalysisResource('orthogroup', payload), false, 'an unbound typed result');
  bindCanonicalTypedResource(orthogroupResult, {
    kind: 'orthogroup-result', type: 'application/json', size: Buffer.byteLength(typedText), lastModified: 0,
    encoding: 'base64', data: Buffer.from(typedText).toString('base64')
  });
  assert.equal(hasRequiredCanonicalAnalysisResource('orthogroup', payload), true, 'a bound typed result');
  assert.equal(hasRequiredCanonicalAnalysisResource('orthogroup', {
    ...payload, pairs: [{ pair_index: 0, tsv: '', rows: [], hit_count: 0 }]
  }), false, 'a pair with only rows');
  // A hit reads the payload back through Vue's reactive cache Map, as a proxy of the bound object.
  const proxied = new Proxy(orthogroupResult, {});
  const previousWindow = globalThis.window;
  globalThis.window = { Vue: { toRaw: (value) => (value === proxied ? orthogroupResult : value) } };
  try {
    assert.equal(hasRequiredCanonicalAnalysisResource('orthogroup', { ...payload, orthogroupResult: proxied }), true);
  } finally {
    globalThis.window = previousWindow;
  }
}
// Current comparison data has one owner; stats can differ without duplicating it.
const immutablePairs = Object.freeze([Object.freeze({ tsv: 'α\tβ\n', canonicalTsv: 'query\tsubject\n' })]);
const immutableResource = Object.freeze({ schema: 3, kind: 'result', value: Object.freeze({
  type: 'CollinearityResult', fields: Object.freeze({ blocks: Object.freeze([]) })
}) });
const conversionPayload = Object.freeze({ pairs: immutablePairs, collinearityResult: immutableResource,
  provenance: Object.freeze({ identity: 'original' }), cache: Object.freeze({ convertedPayloadHit: false }) });
const cachedPayload = stripRuntimeCacheStats(conversionPayload);
assert.deepEqual(cachedPayload, { pairs: immutablePairs, collinearityResult: immutableResource,
  provenance: conversionPayload.provenance });
assert.equal(cachedPayload.pairs, immutablePairs);
assert.equal(cachedPayload.collinearityResult, immutableResource);
assert.equal(conversionPayload.cache.convertedPayloadHit, false);
assert.equal(Object.hasOwn(cachedPayload, 'cache'), false);

// CO-02 (OIPC-C04, PD-OI-022): resolved protein comparisons are reused only
// when the committed Feature visibility rules equal the active rules.
{
  const reuseProbePath = join(identityProbeDir, 'run-analysis-reuse-probe.mjs');
  await writeFile(reuseProbePath, identityProbeSource.replace(
    'const canReuseResolvedProteinArtifacts = ({',
    'export const canReuseResolvedProteinArtifacts = ({'
  ), 'utf8');
  const { canReuseResolvedProteinArtifacts } = await import(pathToFileURL(reuseProbePath));
  const marker = { kind: 'generatedProteinComparison', mode: 'none', settings: {} };
  const comparisons = [marker, { kind: 'orthogroupResult' }];
  const sessionWith = (visibilityText) => ({
    renderRequest: {
      comparisons,
      records: [{ recordKey: 'a', region: {} }],
      diagramOptions: visibilityText === null ? {} : { featureVisibilityTableFile: { resourceId: 'visibility' } }
    },
    resources: visibilityText === null ? {} : {
      visibility: { encoding: 'base64', data: Buffer.from(visibilityText, 'utf8').toString('base64') }
    }
  });
  const active = (featureVisibility) => ({
    featureVisibility, mode: 'orthogroup', candidateLimit: null, bitscore: 50, evalue: 1e-2,
    identity: 0, alignmentLength: 0, orthogroupMembershipMode: 'anchor_core_v1', memberMaxHits: null,
    maxParalogLinks: 2
  });
  const reuse = (session, visibility) => canReuseResolvedProteinArtifacts({
    canonicalComparisons: comparisons, committedSession: session,
    sequences: [{ uid: 'a', region_record_id: '' }], active: active(visibility)
  });
  const hideRule = '*\tCDS\tprotein_id\tCAG38712.1\toff\n';
  assert.equal(reuse(sessionWith(null), ''), true, 'unchanged empty visibility reuses');
  assert.equal(reuse(sessionWith(null), hideRule), false, 'a new hiding rule declines reuse');
  assert.equal(reuse(sessionWith(hideRule), ''), false, 'a removed hiding rule declines reuse');
  assert.equal(reuse(sessionWith(hideRule), hideRule), true, 'the same rules reuse');
  assert.equal(reuse(sessionWith(`${hideRule.trimEnd()}\n\n`), hideRule), true, 'rules compare after normalization');
  assert.equal(reuse({ ...sessionWith(hideRule), resources: {} }, hideRule), false, 'an unreadable committed table declines');

  // OV-1128 B: Similarity groups draws with the paralog links per group, so a changed value declines reuse.
  const reuseParalogs = (committed, current) => canReuseResolvedProteinArtifacts({
    canonicalComparisons: comparisons,
    committedSession: {
      ...sessionWith(null),
      renderRequest: {
        ...sessionWith(null).renderRequest,
        comparisons: [{ ...marker, settings: { collinearMaxParalogLinksPerOrthogroup: committed } }, { kind: 'orthogroupResult' }]
      }
    },
    sequences: [{ uid: 'a', region_record_id: '' }],
    active: { ...active(''), maxParalogLinks: current }
  });
  assert.equal(reuseParalogs(3, 3), true, 'the same paralog links reuse');
  assert.equal(reuseParalogs(3, 2), false, 'changed paralog links decline Similarity groups reuse');
}
