import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';

globalThis.self = {};
if (!globalThis.crypto) globalThis.crypto = webcrypto;

const {
  buildPreparedInputRetention,
  buildPreparedResourceIdentityMap,
  buildGeneratedArtifactTransportIdentity,
  collectGenerationResultTransferList,
  resolveGeneratedFeatureCatalog,
  resolveGenerationCleanupOutcome,
  serializeError
} = await import('../../gbdraw/web/js/workers/diagram-generation-worker.js');
const { readHeldSegments } = await import('../../gbdraw/web/js/services/diagram-worker-protocol.js');

// A resource's cache token is the sha256 hex of its bytes.
const TOKEN_7 = '7'.repeat(64);
const TOKEN_1 = '1'.repeat(64);
const TOKEN_2 = '2'.repeat(64);
const preparedIdentities = buildPreparedResourceIdentityMap([
  {
    resourceId: 'record-1-genbank',
    cacheToken: TOKEN_7,
    size: 1234,
    name: 'secret-name.gbk',
    path: '/temporary/workspace/0001.bin',
    kind: 'genbank'
  }
]);
assert.deepEqual(preparedIdentities, {
  'record-1-genbank': { cacheToken: TOKEN_7, size: 1234 }
});
assert.equal(JSON.stringify(preparedIdentities).includes('secret-name.gbk'), false);
// A4b: Python keeps each drawing's prepared inputs while the retain set names their resources.
assert.deepEqual(buildPreparedInputRetention('linear', [TOKEN_1, TOKEN_2]), {
  drawing: 'linear', resourceTokens: [TOKEN_1, TOKEN_2]
});
assert.deepEqual(buildPreparedInputRetention('circular', undefined), { drawing: 'circular', resourceTokens: [] });
assert.equal(buildPreparedInputRetention(null, [TOKEN_1]), null);
assert.equal(buildPreparedInputRetention(undefined, [TOKEN_1]), null);
for (const drawingId of ['', 'Circular', 'circular/linear', 7]) {
  assert.throws(() => buildPreparedInputRetention(drawingId, []), (error) => error.code === 'HELPER_PROTOCOL');
}
assert.equal(JSON.stringify(preparedIdentities).includes('/temporary/workspace'), false);
for (const cacheToken of ['bad-token', 'render-resource-1']) {
  assert.throws(
    () => buildPreparedResourceIdentityMap([{ resourceId: 'record-1-genbank', cacheToken, size: 1 }]),
    /Invalid cache token/
  );
}

const encoder = new TextEncoder();
const sha256 = async (bytes) => Array.from(
  new Uint8Array(await globalThis.crypto.subtle.digest('SHA-256', bytes)),
  (value) => value.toString(16).padStart(2, '0')
).join('');
const metadataBytes = encoder.encode('{"legendRows":[]}');
const catalogText = '{"schema":4,"items":[]}';
const catalogBytes = encoder.encode(catalogText);

// Step 3b S: a Result is its document pieces; the children texts the page does
// not hold travel in one `segments` buffer.
const CHILDREN = ['<path d="M0,0"/>', '<text>µ</text>'];
const childBytes = CHILDREN.map((text) => encoder.encode(text));
const childRevisions = await Promise.all(childBytes.map(sha256));
const heldList = (indexes) => indexes
  .map((index) => ({ revision: childRevisions[index], bytes: childBytes[index].byteLength }))
  .sort((left, right) => (left.revision < right.revision ? -1 : 1));
const resultOf = (heldIndexes = [], { name = 'out.svg', head = '<svg>' } = {}) => {
  let offset = 0;
  const sent = [];
  const pieces = [{ kind: 'raw', text: '<defs />' }, ...CHILDREN.map((_text, index) => {
    const piece = { kind: 'content', start: `<g id="g${index}">`, revision: childRevisions[index], bytes: childBytes[index].byteLength };
    if (heldIndexes.includes(index)) return piece;
    sent.push(childBytes[index]);
    const placed = { ...piece, offset, length: piece.bytes };
    offset += piece.bytes;
    return placed;
  })];
  const segments = new Uint8Array(offset);
  sent.reduce((at, bytes) => { segments.set(bytes, at); return at + bytes.byteLength; }, 0);
  return { name, document: { head, tail: '</svg>', pieces }, ...(sent.length ? { segments } : {}) };
};
const fullDocument = `<svg><defs />${CHILDREN.map((text, index) => `<g id="g${index}">${text}</g>`).join('')}</svg>`;

// The feature catalog and each Result's segments are their own transferred buffers.
const coldResult = resultOf();
assert.deepEqual(
  collectGenerationResultTransferList({
    results: [coldResult, resultOf([0, 1])],
    metadata: metadataBytes,
    featureCatalog: catalogBytes
  }),
  [coldResult.segments.buffer, metadataBytes.buffer, catalogBytes.buffer]
);
const catalogDigest = { bytes: catalogBytes.byteLength, sha256: await sha256(catalogBytes) };
const identityOf = (overrides = {}, digest = catalogDigest, held = []) => buildGeneratedArtifactTransportIdentity({
  results: [resultOf()],
  metadata: metadataBytes,
  ...overrides
}, digest, held);
const firstIdentity = await identityOf();
const sameIdentity = await identityOf();
const changedIdentity = await identityOf({ results: [resultOf([], { name: 'changed.svg' })] });
const changedContentIdentity = await identityOf({ results: [resultOf([], { head: '<svg id="b">' })] });
const changedMetadataIdentity = await identityOf({ metadata: encoder.encode('{"legendRows":[1]}') });
const changedCatalogIdentity = await identityOf({}, { ...catalogDigest, sha256: '0'.repeat(64) });
const noCatalogIdentity = await identityOf({}, null);
// Held or sent, the same document is the same artifact with the same byte counts (B3).
const partlyHeldIdentity = await identityOf({ results: [resultOf([1])] }, catalogDigest, heldList([1]));
const allHeldIdentity = await identityOf({ results: [resultOf([0, 1])] }, catalogDigest, heldList([0, 1]));
assert.match(firstIdentity.fingerprint, /^[0-9a-f]{64}$/);
assert.equal(firstIdentity.fingerprint, sameIdentity.fingerprint);
assert.deepEqual(partlyHeldIdentity, firstIdentity);
assert.deepEqual(allHeldIdentity, firstIdentity);
assert.notEqual(firstIdentity.fingerprint, changedIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, changedContentIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, changedMetadataIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, changedCatalogIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, noCatalogIdentity.fingerprint);
assert.equal(firstIdentity.schema, 1);
assert.equal(firstIdentity.resultBytes, encoder.encode(fullDocument).byteLength);
assert.equal(firstIdentity.metadataBytes, metadataBytes.byteLength + catalogBytes.byteLength);
assert.equal(
  firstIdentity.retainedBytes,
  encoder.encode(fullDocument).byteLength
    + (metadataBytes.byteLength + catalogBytes.byteLength) * 2
    + encoder.encode('out.svg').byteLength
);

// A reply that does not match its request or its own bytes is a classified
// result-admission error: a held text the request did not name, a sent text
// whose sha256 or length differs, an unknown piece kind, a gap in the buffer.
const rejectsReply = async (result, held = []) => {
  await assert.rejects(identityOf({ results: [result] }, catalogDigest, held), (error) => {
    assert.equal(error.code, 'RESULT_INVALID');
    assert.equal(error.stage, 'result-admission');
    return true;
  });
};
await rejectsReply(resultOf([1]));
await rejectsReply(resultOf([1]), heldList([0]));
await rejectsReply(resultOf([1]), [{ ...heldList([1])[0], bytes: 1 }]);
const badSha = resultOf();
badSha.document.pieces[1] = { ...badSha.document.pieces[1], revision: '0'.repeat(64) };
await rejectsReply(badSha);
const badLength = resultOf();
badLength.document.pieces[2] = { ...badLength.document.pieces[2], bytes: 3, length: 3 };
await rejectsReply(badLength);
const unknownKind = resultOf();
unknownKind.document.pieces[0] = { kind: 'comment', text: '<!-- x -->' };
await rejectsReply(unknownKind);
const withGap = resultOf();
withGap.segments = new Uint8Array([...withGap.segments, 32]);
await rejectsReply(withGap);
await rejectsReply({ name: 'out.svg', content: encoder.encode('<svg />') });

// The run message's held segments: sorted, distinct, sha256 revisions with byte counts.
const invalidHeld = () => Object.assign(new Error('invalid'), { code: 'HELPER_PROTOCOL' });
assert.deepEqual(readHeldSegments(null, invalidHeld), []);
assert.deepEqual(readHeldSegments(heldList([0, 1]), invalidHeld), heldList([0, 1]));
for (const malformed of [
  {}, [{ revision: 'e:1', bytes: 1 }], [{ revision: childRevisions[0], bytes: -1 }],
  [...heldList([0, 1])].reverse(), [...heldList([0]), ...heldList([0])]
]) assert.throws(() => readHeldSegments(malformed, invalidHeld), /invalid/);

// A sent catalog is hashed; "unchanged" takes the held digest, so both name
// one artifact; without a catalog and a held digest there is no revision.
const sent = await resolveGeneratedFeatureCatalog({
  results: [resultOf()],
  metadata: metadataBytes,
  featureCatalog: catalogBytes,
  featureCatalogRevision: 'e:1'
}, null, []);
assert.deepEqual(sent.featureCatalogRevision, { revision: 'e:1', ...catalogDigest });
assert.equal(sent.artifactIdentity.fingerprint, firstIdentity.fingerprint);
const unchanged = await resolveGeneratedFeatureCatalog({
  results: [resultOf([0])],
  metadata: metadataBytes,
  featureCatalog: undefined,
  featureCatalogRevision: 'e:1'
}, { revision: 'e:1', ...catalogDigest }, heldList([0]));
assert.deepEqual(unchanged.featureCatalogRevision, { revision: 'e:1', ...catalogDigest });
assert.equal(unchanged.artifactIdentity.fingerprint, firstIdentity.fingerprint);
const unnamed = await resolveGeneratedFeatureCatalog({
  results: [resultOf()],
  metadata: metadataBytes,
  featureCatalog: catalogBytes,
  featureCatalogRevision: undefined
}, null, []);
assert.equal(unnamed.featureCatalogRevision, null);
assert.equal(unnamed.artifactIdentity.fingerprint, firstIdentity.fingerprint);
const neither = await resolveGeneratedFeatureCatalog({
  results: [resultOf()],
  metadata: metadataBytes,
  featureCatalogRevision: 'e:1'
}, null, []);
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

// OV-1021: only a failed snapshot step asks the page for a Worker that boots
// cold (snapshotError); any other init error is an ordinary init error, also
// when the init message offers a snapshot.
{
  /** @param {Record<string, unknown>} message */
  const initReply = (message) => new Promise((resolve) => {
    self.postMessage = (reply) => { if (reply.type === 'init') resolve(reply); };
    self.onmessage({ data: { type: 'init', ...message } });
  });
  const offline = 'http://127.0.0.1:9';
  const init = {
    pyodideIndexUrl: `${offline}/pyodide/`,
    gbdrawWheelUrl: `${offline}/gbdraw.whl`,
    pyodideSnapshotUrl: `${offline}/snapshot/`
  };
  const originalFetch = globalThis.fetch;
  try {
    // The snapshot download is still pending when the Pyodide module fails.
    globalThis.fetch = () => new Promise(() => {});
    const brokenModule = await initReply({
      ...init, id: 'broken-module', pyodideModuleUrl: 'data:text/javascript,throw new Error("broken Pyodide module")'
    });
    assert.equal(brokenModule.ok, false);
    assert.equal('snapshotError' in brokenModule, false);
    assert.equal(brokenModule.error.code, 'WORKER_INIT');

    globalThis.fetch = () => Promise.reject(new TypeError('fetch failed'));
    const unreachableSnapshot = await initReply({
      ...init, id: 'unreachable-snapshot', pyodideModuleUrl: 'data:text/javascript,export const loadPyodide = () => null;'
    });
    assert.equal(unreachableSnapshot.ok, false);
    assert.equal(unreachableSnapshot.snapshotError, 'fetch failed');
  } finally {
    globalThis.fetch = originalFetch;
  }
}

// A cold boot starts every wheel download (one GET each, no HEAD check)
// before Pyodide loads.
{
  const events = [];
  const offline = 'http://127.0.0.1:9';
  const localWheelUrls = [`${offline}/a-1.0-py3-none-any.whl`, `${offline}/b-1.0-py3-none-any.whl`];
  const gbdrawWheelUrl = `${offline}/gbdraw-0.0.0-py3-none-any.whl?v=1`;
  globalThis.__FAKE_LOAD_PYODIDE__ = async () => {
    events.push(['loadPyodide']);
    return { loadPackage: async () => { throw new Error('stop after loadPyodide'); } };
  };
  const originalFetch = globalThis.fetch;
  globalThis.fetch = async (url, options = {}) => {
    events.push([options.method || 'GET', String(url)]);
    return new Response('wheel');
  };
  try {
    const reply = await new Promise((resolve) => {
      self.postMessage = (message) => { if (message.type === 'init') resolve(message); };
      self.onmessage({ data: {
        type: 'init', id: 'cold-boot-downloads', pyodideIndexUrl: `${offline}/pyodide/`, localWheelUrls, gbdrawWheelUrl,
        pyodideModuleUrl: 'data:text/javascript,export const loadPyodide = (options) => globalThis.__FAKE_LOAD_PYODIDE__(options);'
      } });
    });
    assert.equal(reply.ok, false);
    assert.deepEqual(events, [...[...localWheelUrls, gbdrawWheelUrl].map((url) => ['GET', url]), ['loadPyodide']]);
  } finally {
    globalThis.fetch = originalFetch;
    delete globalThis.__FAKE_LOAD_PYODIDE__;
  }
}

// Every helper that stages resources accepts each key the page's resource
// staging writes into its payload (PR #1313 CI: `extractCdsProteinFasta` and
// `convertStoredProteinResultFrames` refused `retainedResourceTokens`).
{
  const { helperPayloadKeys } = await import('../../gbdraw/web/js/workers/diagram-generation-worker.js');
  const { DIAGRAM_HELPER_OPERATIONS, DIAGRAM_HELPER_STAGING_OPERATIONS } = await import(
    '../../gbdraw/web/js/services/diagram-worker-protocol.js'
  );
  const { createDiagramResourceTransport } = await import('../../gbdraw/web/js/services/diagram-resource-staging.js');
  const prepared = await createDiagramResourceTransport().prepare({
    request: { records: [], diagramOptions: {} }, resources: {}, retain: true
  });
  const stagedKeys = Object.keys(prepared).filter((key) => typeof prepared[key] !== 'function');
  assert.ok(stagedKeys.includes('retainedResourceTokens'));
  assert.ok(DIAGRAM_HELPER_STAGING_OPERATIONS.includes(DIAGRAM_HELPER_OPERATIONS.EXTRACT_CDS_PROTEIN_FASTA));
  for (const operation of DIAGRAM_HELPER_STAGING_OPERATIONS) {
    const accepted = helperPayloadKeys(operation);
    for (const key of stagedKeys) assert.ok(accepted?.includes(key), `${operation} accepts ${key}`);
  }
  assert.equal(helperPayloadKeys(DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES)?.includes('stagedResources'), false);
  assert.equal(helperPayloadKeys('notAHelper'), null);
}
