import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';

globalThis.self = {};
if (!globalThis.crypto) globalThis.crypto = webcrypto;

const {
  buildPreparedInputRetention,
  buildPreparedResourceIdentityMap,
  buildGeneratedArtifactTransportIdentity,
  callJsonHelper,
  collectGenerationResultTransferList,
  createReplyBuffer,
  resolveGeneratedFeatureCatalog,
  nameStagedHelperInput,
  resolveGenerationCleanupOutcome,
  serializeError,
  noteDroppedReply,
  cancelDrainReply
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
assert.deepEqual(buildPreparedInputRetention('linear', [TOKEN_1, TOKEN_2], []), {
  drawing: 'linear', resourceTokens: [TOKEN_1, TOKEN_2], heldTokens: []
});
assert.deepEqual(buildPreparedInputRetention('circular', undefined, []), { drawing: 'circular', resourceTokens: [], heldTokens: [] });
// A render message without `heldResourceTokens` is a protocol error, not an empty list.
assert.throws(
  () => buildPreparedInputRetention('circular', [TOKEN_1]),
  (error) => error.code === 'HELPER_PROTOCOL' && error.context.reason === 'retention-held'
);
assert.equal(buildPreparedInputRetention(null, [TOKEN_1]), null);
assert.equal(buildPreparedInputRetention(undefined, [TOKEN_1]), null);
for (const drawingId of ['', 'Circular', 'circular/linear', 7]) {
  assert.throws(() => buildPreparedInputRetention(drawingId, [], []), (error) => error.code === 'HELPER_PROTOCOL');
}
// Held tokens reach the prepared input retention; a held token the Worker does not keep is a protocol error.
assert.deepEqual(buildPreparedInputRetention('circular', [TOKEN_1, TOKEN_2], [TOKEN_2]), {
  drawing: 'circular', resourceTokens: [TOKEN_1, TOKEN_2], heldTokens: [TOKEN_2]
});
assert.throws(
  () => buildPreparedInputRetention('circular', [TOKEN_1], [TOKEN_2]),
  (error) => error.code === 'HELPER_PROTOCOL' && error.context.reason === 'retention-held'
);
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

// Python writes each chunk of a reply buffer to the buffer's file; `take` moves
// it out of the file system, and the buffer joins the chunks once.
const replyFiles = () => {
  const files = new Map();
  return {
    files,
    readFile: (path) => new Uint8Array(files.get(path)),
    unlink: (path) => { assert.ok(files.delete(path)); }
  };
};
{
  const fs = replyFiles();
  const buffer = createReplyBuffer(fs);
  for (const chunk of [[1, 2, 3], [], [4, 5, 6, 7]]) {
    fs.files.set(buffer.path, chunk);
    buffer.take();
    assert.equal(fs.files.size, 0);
  }
  const joined = buffer.finish();
  assert.ok(joined instanceof ArrayBuffer);
  assert.deepEqual([...new Uint8Array(joined)], [1, 2, 3, 4, 5, 6, 7]);
  assert.equal(createReplyBuffer(fs).finish().byteLength, 0);
}

// A sequence read replies in UTF-8 through a reply buffer (`_replied_in_utf8` in python-helpers.js): the helper
// returns the buffer, and its bytes are the reply's JSON.
{
  const fs = replyFiles();
  const text = JSON.stringify({ spans: [{ sequence: 'ACGT', recordId: 'é\u2028' }] });
  const helper = (name, newReplyBuffer, ...rest) => {
    assert.equal(name, 'read_match_sequences_json');
    assert.deepEqual(rest, ['{}']);
    const buffer = newReplyBuffer();
    fs.files.set(buffer.path, [...new TextEncoder().encode(text)]);
    buffer.take();
    return buffer;
  };
  helper.destroy = () => {};
  assert.deepEqual(
    callJsonHelper({ globals: { get: () => helper } }, 'read_match_sequences_json', [() => createReplyBuffer(fs), '{}']),
    JSON.parse(text)
  );
}

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
    + metadataBytes.byteLength + catalogBytes.byteLength
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
  // Only a Generate sends `heldResourceTokens`: a helper keeps every file for the next Generate.
  const stagedKeys = Object.keys(prepared).filter((key) => typeof prepared[key] !== 'function' && key !== 'heldResourceTokens');
  assert.ok(stagedKeys.includes('retainedResourceTokens'));
  assert.ok(DIAGRAM_HELPER_STAGING_OPERATIONS.includes(DIAGRAM_HELPER_OPERATIONS.EXTRACT_CDS_PROTEIN_FASTA));
  for (const operation of DIAGRAM_HELPER_STAGING_OPERATIONS) {
    const accepted = helperPayloadKeys(operation);
    for (const key of stagedKeys) assert.ok(accepted?.includes(key), `${operation} accepts ${key}`);
  }
  assert.equal(helperPayloadKeys(DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES)?.includes('stagedResources'), false);
  assert.equal(helperPayloadKeys('notAHelper'), null);
  // Load of a Session 27-44: the one Worker call takes the document, whose
  // resources travel inside it, and the schema of the request the file
  // stored (OV-2052), and stages no file (DESIGN 4.3).
  assert.deepEqual(helperPayloadKeys(DIAGRAM_HELPER_OPERATIONS.UPGRADE_LEGACY_SESSION), ['document', 'storedRequestSchema']);
  assert.equal(DIAGRAM_HELPER_STAGING_OPERATIONS.includes(DIAGRAM_HELPER_OPERATIONS.UPGRADE_LEGACY_SESSION), false);
  // New drawing from a region: Python derives the drawing from the source
  // drawing's draft projection, whose resources the page stages.
  const { stagedHelperReferences } = await import('../../gbdraw/web/js/services/diagram-worker-protocol.js');
  const region = DIAGRAM_HELPER_OPERATIONS.DERIVE_REGION_DRAWING;
  assert.equal(region, 'deriveRegionDrawing');
  assert.ok(DIAGRAM_HELPER_STAGING_OPERATIONS.includes(region));
  for (const key of ['canonicalRequest', 'regions', 'adaptSizes']) assert.ok(helperPayloadKeys(region)?.includes(key), key);
  const canonicalRequest = { records: [] };
  assert.equal(stagedHelperReferences(region, { canonicalRequest }), canonicalRequest);
}

// D7: a Cancel drains the queue. The interrupt and then Python's cancel flag
// are cleared after the operation in flight ends, the replies of every request
// up to `through` are dropped (a pending part acknowledgement is released),
// and the answer names the Cancel.
{
  const offline = 'http://127.0.0.1:9';
  const replies = [];
  self.postMessage = (message) => { replies.push(message); };
  const settle = async () => { for (let i = 0; i < 50; i += 1) await new Promise((resolve) => setImmediate(resolve)); };
  // A run whose reply waits on work outside the event loop (the joined reply buffers) takes a varying
  // number of turns: wait for its run reply, bounded, instead of a fixed count.
  const settleUntilRun = async () => {
    for (let i = 0; i < 2000 && !replies.some(({ type }) => type === 'run'); i += 1) {
      await new Promise((resolve) => (i < 50 ? setImmediate(resolve) : setTimeout(resolve, 1)));
    }
  };
  const bootMessage = (id, interruptBuffer) => ({
    type: 'init', id, interruptBuffer, pyodideIndexUrl: `${offline}/pyodide/`, localWheelUrls: [],
    gbdrawWheelUrl: `${offline}/gbdraw-0.0.0-py3-none-any.whl`,
    pyodideModuleUrl: 'data:text/javascript,export const loadPyodide = (options) => globalThis.__FAKE_LOAD_PYODIDE__(options);'
  });
  const originalFetch = globalThis.fetch;
  globalThis.fetch = async () => new Response('wheel');
  try {
    // Without a runtime (its boot is in flight, then fails): the buffer is cleared, after the boot ends.
    const buffer = new Uint8Array(new SharedArrayBuffer(1));
    let failBoot;
    globalThis.__FAKE_LOAD_PYODIDE__ = () => new Promise((_resolve, reject) => { failBoot = reject; });
    self.onmessage({ data: bootMessage('boot-in-flight', buffer) });
    buffer[0] = 2;
    self.onmessage({ data: { type: 'cancel', id: 'cancel-during-boot', through: 0 } });
    await settle();
    assert.deepEqual(replies.map(({ type }) => type), ['init-progress']);
    failBoot(new Error('stop the boot'));
    await settle();
    assert.deepEqual(replies.filter(({ type }) => type !== 'init-progress').map(({ type, ok }) => [type, ok]),
      [['init', false], ['cancel', true]]);
    assert.equal(replies.at(-1).id, 'cancel-during-boot');
    assert.equal(buffer[0], 0);

    // A booted runtime: the buffer reaches Pyodide once, at boot.
    replies.length = 0;
    const helperCalls = [];
    let interruptBuffer = null;
    let helperReply = () => 'null';
    let wrapperNotice = null;
    let wrapperReply = null;
    // The file each reply chunk crosses through (ReplyBuffer.take).
    const chunkFiles = new Map();
    const fs = new Proxy({
      analyzePath: () => ({ exists: false }),
      readdir: () => [],
      readFile: (path) => chunkFiles.get(path),
      unlink: (path) => { chunkFiles.delete(path); }
    }, {
      get: (target, name) => target[name] ?? (() => {})
    });
    const pyodide = {
      loadPackage: async () => {},
      pyimport: () => ({ install: async () => {} }),
      FS: fs,
      runPythonAsync: async () => {},
      runPython: (code) => (code.includes('get_web_runtime_capabilities') ? '{}' : undefined),
      setInterruptBuffer: (value) => { interruptBuffer = value; },
      globals: {
        get: (name) => ({
          call_web_json_helper: (name, ...args) => {
            helperCalls.push(name);
            return helperReply(name, ...args);
          },
          // CAPACITY-WARN: Python calls its notice argument (the one before the reply buffer
          // factory) with each sentence while it renders.
          run_canonical_request_wrapper: (...args) => {
            assert.equal(typeof args.at(-1), 'function');
            if (wrapperNotice) args.at(-2)(wrapperNotice);
            if (wrapperReply) return wrapperReply(args.at(-1));
            return { error: { code: 'INPUT_INVALID', operation: 'generate', stage: 'render' } };
          }
        })[name]
      }
    };
    globalThis.__FAKE_LOAD_PYODIDE__ = async () => pyodide;
    const booted = new Uint8Array(new SharedArrayBuffer(1));
    self.onmessage({ data: bootMessage('booted', booted) });
    await settle();
    assert.equal(replies.find(({ type }) => type === 'init')?.ok, true);
    assert.equal(interruptBuffer, booted);

    // A run that is not canceled posts its result; a canceled one posts nothing.
    replies.length = 0;
    const run = (requestId) => self.onmessage({ data: { type: 'run', requestId, payload: {
      request: {}, resourceManifest: [], stagedResources: [], retainedResourceTokens: []
    } } });
    run(10);
    await settle();
    assert.deepEqual(replies.filter(({ type }) => type === 'run').map(({ requestId, ok }) => [requestId, ok]), [[10, true]]);
    // A notice Python sends while it renders goes out as a progress message of that run, before its result.
    replies.length = 0;
    wrapperNotice = JSON.stringify({ recordId: 'NC_1', labelKind: 'circular_outer', labelCount: 195, fontSize: 14, capacityRatio: 1.4 });
    run(13);
    await settle();
    assert.deepEqual(replies.filter(({ type }) => type !== 'test-lifecycle').map(({ type, requestId, stage, notice }) => (
      [type, requestId, stage, notice]
    )), [
      ['progress', 13, 'preparing-resources', undefined],
      ['progress', 13, 'rendering', undefined],
      ['progress', 13, 'rendering', wrapperNotice],
      ['progress', 13, 'finalizing', undefined],
      ['run', 13, undefined, undefined]
    ]);
    wrapperNotice = null;
    // The reply buffers Python wrote are joined after it returns, at result admission:
    // the run carries exact-length buffers, and a buffer too large to allocate fails
    // there, as the conversion of Python's bytes did.
    const writeChunk = (buffer, text) => {
      chunkFiles.set(buffer.path, new TextEncoder().encode(text));
      buffer.take();
    };
    replies.length = 0;
    wrapperReply = (newReplyBuffer) => {
      const metadata = newReplyBuffer();
      writeChunk(metadata, '{"a":');
      writeChunk(metadata, '1}');
      return { results: [], metadata };
    };
    run(14);
    await settleUntilRun();
    const joinedRun = replies.find(({ type }) => type === 'run');
    assert.equal(joinedRun.ok, true);
    assert.equal(new TextDecoder().decode(joinedRun.results.metadata), '{"a":1}');
    replies.length = 0;
    wrapperReply = (newReplyBuffer) => {
      const metadata = newReplyBuffer();
      metadata.finish = () => { throw new RangeError('Array buffer allocation failed'); };
      return { results: [], metadata };
    };
    run(15);
    await settleUntilRun();
    const failedRun = replies.find(({ type }) => type === 'run');
    assert.deepEqual([failedRun.ok, failedRun.error.code, failedRun.error.stage], [false, 'UNKNOWN', 'result-admission']);
    wrapperReply = null;
    replies.length = 0;
    run(11);
    booted[0] = 2;
    self.onmessage({ data: { type: 'cancel', id: 'cancel-run', through: 11 } });
    await settle();
    // Progress may still go out (the page ignores it); no result does.
    assert.deepEqual(replies.filter(({ type }) => type !== 'progress').map(({ type, id }) => [type, id]),
      [['cancel', 'cancel-run']]);
    assert.equal(booted[0], 0);
    assert.equal(helperCalls.at(-1), 'clear_web_cancel');

    // A bounded reply waiting for its part acknowledgement: the Cancel releases it, and the rest is dropped.
    replies.length = 0;
    helperReply = () => JSON.stringify({ rows: Array.from({ length: 3000 }, (_, index) => `${index}`.padEnd(100, 'x')) });
    self.onmessage({ data: { type: 'helper', requestId: 12, operation: 'evaluateRules', payload: { rules: [] } } });
    await settle();
    assert.deepEqual(replies.map(({ type, status }) => [type, status]), [['helper', 'part']]);
    self.onmessage({ data: { type: 'cancel', id: 'cancel-helper', through: 12 } });
    await settle();
    assert.deepEqual(replies.map(({ type, status, id }) => [type, status, id]),
      [['helper', 'part', undefined], ['cancel', undefined, 'cancel-helper']]);
  } finally {
    globalThis.fetch = originalFetch;
    delete globalThis.__FAKE_LOAD_PYODIDE__;
  }
}

// OV-3027: a canceled request's reply is dropped, but when it says the engine ran out of memory the drain answer
// carries that error once, so the page still replaces the engine.
{
  const outOfMemory = { code: 'OUT_OF_MEMORY', operation: 'generate', stage: 'render', context: {} };
  assert.deepEqual(cancelDrainReply('a'), { type: 'cancel', id: 'a', ok: true });
  noteDroppedReply({ code: 'RENDER_FAILED', operation: 'generate', stage: 'render', context: {} });
  noteDroppedReply(undefined);
  assert.deepEqual(cancelDrainReply('b'), { type: 'cancel', id: 'b', ok: true });
  noteDroppedReply(outOfMemory);
  assert.deepEqual(cancelDrainReply('c'), { type: 'cancel', id: 'c', ok: true, error: outOfMemory });
  assert.deepEqual(cancelDrainReply('d'), { type: 'cancel', id: 'd', ok: true });
}

// PD-OI-116: a helper reads a chosen gzip file as chosen, and Python refuses one over the cap by the
// name of the file it read, which the Worker staged under its role's name; the refusal names the file
// the user chose instead, with the cap. A name that is not a staged one stays.
{
  const refusal = serializeError({ code: 'INPUT_UNREADABLE', operation: 'unknown', stage: 'helper',
    context: { reason: 'GZIP_EXPANDED_SIZE', limit: 512, inputName: 'source.bin' } });
  const named = nameStagedHelperInput(refusal, [{ role: 'source', name: 'big genome.fna.gz', blob: new Blob(['x']) }]);
  assert.equal(named.context.inputName, 'big genome.fna.gz');
  assert.match(named.summary, /big genome\.fna\.gz/);
  assert.match(named.summary, /512 MiB/);
  assert.doesNotMatch(named.summary, /source\.bin/);
  const unnamed = nameStagedHelperInput(refusal, [{ role: 'source', blob: new Blob(['x']) }]);
  assert.equal(unnamed.context.inputName, undefined, 'a staged name is never shown');
  const generate = serializeError({ code: 'INPUT_UNREADABLE', stage: 'render',
    context: { reason: 'GZIP_EXPANDED_SIZE', limit: 512, inputName: 'chosen.gb.gz' } });
  assert.equal(nameStagedHelperInput(generate, []), generate);
}
