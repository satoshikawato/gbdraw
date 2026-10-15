import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';

globalThis.self = {};
if (!globalThis.crypto) globalThis.crypto = webcrypto;

const {
  buildPreparedResourceIdentityMap,
  buildGeneratedArtifactTransportIdentity,
  collectGenerationResultTransferList,
  resolveGeneratedFeatureCatalog,
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
const metadataBytes = new TextEncoder().encode('{"legendRows":[]}');
const catalogText = '{"schema":4,"items":[]}';
const catalogBytes = new TextEncoder().encode(catalogText);
// Step 3a: the feature catalog is its own transferred buffer beside the metadata.
assert.deepEqual(
  collectGenerationResultTransferList({
    results: [{ content: svgBytes }, { content: '<svg />' }],
    metadata: metadataBytes,
    featureCatalog: catalogBytes
  }),
  [svgBytes.buffer, metadataBytes.buffer, catalogBytes.buffer]
);
const catalogDigest = {
  bytes: catalogBytes.byteLength,
  sha256: Array.from(
    new Uint8Array(await globalThis.crypto.subtle.digest('SHA-256', catalogBytes)),
    (value) => value.toString(16).padStart(2, '0')
  ).join('')
};
const identityOf = (overrides = {}, digest = catalogDigest) => buildGeneratedArtifactTransportIdentity({
  results: [{ name: 'out.svg', content: svgBytes }],
  metadata: metadataBytes,
  ...overrides
}, digest);
const firstIdentity = await identityOf();
const sameIdentity = await identityOf();
const changedIdentity = await identityOf({ results: [{ name: 'changed.svg', content: svgBytes }] });
const changedContentIdentity = await identityOf({
  results: [{ name: 'out.svg', content: new TextEncoder().encode('<svg>b</svg>') }]
});
const changedMetadataIdentity = await identityOf({ metadata: new TextEncoder().encode('{"legendRows":[1]}') });
const changedCatalogIdentity = await identityOf({}, { ...catalogDigest, sha256: '0'.repeat(64) });
const noCatalogIdentity = await identityOf({}, null);
assert.match(firstIdentity.fingerprint, /^[0-9a-f]{64}$/);
assert.equal(firstIdentity.fingerprint, sameIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, changedIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, changedContentIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, changedMetadataIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, changedCatalogIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, noCatalogIdentity.fingerprint);
assert.equal(firstIdentity.schema, 1);
assert.equal(firstIdentity.metadataBytes, metadataBytes.byteLength + catalogBytes.byteLength);
assert.equal(
  firstIdentity.retainedBytes,
  svgBytes.byteLength
    + (metadataBytes.byteLength + catalogBytes.byteLength) * 2
    + new TextEncoder().encode('out.svg').byteLength
);

// A sent catalog is hashed; "unchanged" takes the held digest, so both name
// one artifact; without a catalog and a held digest there is no revision.
const sent = await resolveGeneratedFeatureCatalog({
  results: [{ name: 'out.svg', content: svgBytes }],
  metadata: metadataBytes,
  featureCatalog: catalogBytes,
  featureCatalogRevision: 'e:1'
}, null);
assert.deepEqual(sent.featureCatalogRevision, { revision: 'e:1', ...catalogDigest });
assert.equal(sent.artifactIdentity.fingerprint, firstIdentity.fingerprint);
const unchanged = await resolveGeneratedFeatureCatalog({
  results: [{ name: 'out.svg', content: svgBytes }],
  metadata: metadataBytes,
  featureCatalog: undefined,
  featureCatalogRevision: 'e:1'
}, { revision: 'e:1', ...catalogDigest });
assert.deepEqual(unchanged.featureCatalogRevision, { revision: 'e:1', ...catalogDigest });
assert.equal(unchanged.artifactIdentity.fingerprint, firstIdentity.fingerprint);
const unnamed = await resolveGeneratedFeatureCatalog({
  results: [{ name: 'out.svg', content: svgBytes }],
  metadata: metadataBytes,
  featureCatalog: catalogBytes,
  featureCatalogRevision: undefined
}, null);
assert.equal(unnamed.featureCatalogRevision, null);
assert.equal(unnamed.artifactIdentity.fingerprint, firstIdentity.fingerprint);
const neither = await resolveGeneratedFeatureCatalog({
  results: [{ name: 'out.svg', content: svgBytes }],
  metadata: metadataBytes,
  featureCatalogRevision: 'e:1'
}, null);
assert.equal(neither.featureCatalogRevision, null);
assert.equal(neither.artifactIdentity.fingerprint, noCatalogIdentity.fingerprint);

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
