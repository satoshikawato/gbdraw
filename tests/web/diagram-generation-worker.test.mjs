import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';

globalThis.self = {};
if (!globalThis.crypto) globalThis.crypto = webcrypto;

const {
  buildPreparedResourceIdentityMap,
  buildGeneratedArtifactTransportIdentity,
  collectGenerationResultTransferList,
  resolveGenerationCleanupOutcome,
  serializeError
} = await import('../../gbdraw/web/js/workers/diagram-generation-worker.js');

const preparedIdentities = buildPreparedResourceIdentityMap([
  {
    resourceId: 'record-1-genbank',
    cacheToken: 'render-resource-7',
    size: 1234,
    name: 'secret-name.gbk',
    path: '/temporary/workspace/0001.bin',
    kind: 'genbank'
  }
]);
assert.deepEqual(preparedIdentities, {
  'record-1-genbank': { cacheToken: 'render-resource-7', size: 1234 }
});
assert.equal(JSON.stringify(preparedIdentities).includes('secret-name.gbk'), false);
assert.equal(JSON.stringify(preparedIdentities).includes('/temporary/workspace'), false);
assert.throws(
  () => buildPreparedResourceIdentityMap([
    { resourceId: 'record-1-genbank', cacheToken: 'bad-token', size: 1 }
  ]),
  /Invalid cache token/
);
assert.throws(
  () => buildPreparedResourceIdentityMap([
    { resourceId: 'record-1-genbank', cacheToken: 'render-resource-1', size: 1 },
    { resourceId: 'record-1-genbank', cacheToken: 'render-resource-2', size: 1 }
  ]),
  /duplicate/
);

const svgBytes = new TextEncoder().encode('<svg />');
const metadataBytes = new TextEncoder().encode('{"featureCatalog":{}}');
assert.deepEqual(
  collectGenerationResultTransferList({
    results: [{ content: svgBytes }, { content: '<svg />' }],
    metadata: metadataBytes
  }),
  [svgBytes.buffer, metadataBytes.buffer]
);
const firstIdentity = await buildGeneratedArtifactTransportIdentity({
  results: [{ name: 'out.svg', content: svgBytes }],
  metadata: metadataBytes
});
const sameIdentity = await buildGeneratedArtifactTransportIdentity({
  results: [{ name: 'out.svg', content: svgBytes }],
  metadata: metadataBytes
});
const changedIdentity = await buildGeneratedArtifactTransportIdentity({
  results: [{ name: 'changed.svg', content: svgBytes }],
  metadata: metadataBytes
});
const changedContentIdentity = await buildGeneratedArtifactTransportIdentity({
  results: [{ name: 'out.svg', content: new TextEncoder().encode('<svg>b</svg>') }],
  metadata: metadataBytes
});
const changedMetadataIdentity = await buildGeneratedArtifactTransportIdentity({
  results: [{ name: 'out.svg', content: svgBytes }],
  metadata: new TextEncoder().encode('{"featureCatalog":{"schema":4}}')
});
assert.match(firstIdentity.fingerprint, /^[0-9a-f]{64}$/);
assert.equal(firstIdentity.fingerprint, sameIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, changedIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, changedContentIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, changedMetadataIdentity.fingerprint);
assert.equal(
  firstIdentity.retainedBytes,
  svgBytes.byteLength
    + metadataBytes.byteLength * 2
    + new TextEncoder().encode('out.svg').byteLength
);

const sentinel = 'PRIVATE_SENTINEL';
const pythonPayload = { error: { code: 'COMPARISON_IDENTITY', operation: 'generate',
  stage: 'render', context: { reason: 'SOURCE_VIEW_CONFLICT' } } };
const pythonResult = resolveGenerationCleanupOutcome({
  result: pythonPayload, destroyError: new Error(sentinel), workspaceError: new Error(sentinel)
});
assert.equal(pythonResult, pythonPayload);
assert.deepEqual(pythonResult.error.secondary, [
  { code: 'CLEANUP_FAILED', stage: 'cleanup' }, { code: 'CLEANUP_FAILED', stage: 'cleanup' }
]);
assert.equal(serializeError(pythonResult.error).code, 'COMPARISON_IDENTITY');
assert.doesNotMatch(JSON.stringify(serializeError(pythonResult.error)), /PRIVATE_SENTINEL/);

const primaryJsError = Object.assign(new Error(sentinel), { code: 'RESOURCE_INVALID', stage: 'resource-staging' });
assert.throws(() => resolveGenerationCleanupOutcome({ primaryError: primaryJsError,
  destroyError: new Error(sentinel), workspaceError: new Error(sentinel) }), (error) => {
  assert.equal(error, primaryJsError);
  const serialized = serializeError(error);
  assert.equal(serialized.code, 'RESOURCE_INVALID');
  assert.equal(serialized.secondary.length, 2);
  assert.doesNotMatch(JSON.stringify(serialized), /PRIVATE_SENTINEL/);
  return true;
});
const destroyError = new Error(sentinel);
assert.throws(() => resolveGenerationCleanupOutcome({ result: { results: [] }, destroyError,
  workspaceError: new Error(sentinel) }), (error) => {
  assert.equal(error, destroyError);
  assert.equal(serializeError(error).code, 'CLEANUP_FAILED');
  assert.equal(serializeError(error).stage, 'cleanup');
  assert.equal(serializeError(error).secondary.length, 1);
  return true;
});
const workspaceError = new Error(sentinel);
assert.throws(() => resolveGenerationCleanupOutcome({ result: { results: [] }, workspaceError }),
  (error) => error === workspaceError && serializeError(error).code === 'CLEANUP_FAILED');
