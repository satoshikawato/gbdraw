import assert from 'node:assert/strict';
import { existsSync, readFileSync } from 'node:fs';

const {
  MAX_EXPANDED_SESSION_BYTES,
  largeSessionDownloadMegabytes,
  sessionContainerResources,
  writeSessionFile
} = await import('../../gbdraw/web/js/services/session-file.js');
const { openSessionContainerEntry, readSessionContainer } = await import('../../gbdraw/web/js/services/session-container.js');
const { askBeforeLargeDownload, downloadBlob } = await import('../../gbdraw/web/js/services/text-download.js');

// A Save that keeps every entry in memory: the browser lets no copy be written.
const NO_COPIES = { keep: null };
const bytesOf = (text) => new TextEncoder().encode(text);
const streamOf = (bytes) => new Blob([bytes]).stream();
/** The Session file of a document with no inputs. */
const writeDocument = (value, sources = NO_COPIES) => writeSessionFile(value, {
  readResource: async () => assert.fail('a document with no inputs reads none'), ...sources
});
/** The text of `gbdraw-session.json` in the file Save writes for `value`. */
const documentText = async (value) => {
  const { blob } = await writeDocument(value);
  const { documentEntry } = await readSessionContainer(blob);
  return new Response(openSessionContainerEntry(blob, documentEntry)).text();
};
const assertJsonParity = async (value) => {
  const actual = await documentText(value);
  assert.equal(actual, JSON.stringify(value));
  return JSON.parse(actual);
};

// Save writes one Session file container: `gbdraw-session.json`, then one entry per
// input, read from its source by resource id, gzip bytes stored and others deflated.
const genome = bytesOf('LOCUS       test\n'.repeat(500));
const gzipInput = new Uint8Array([0x1f, 0x8b, 8, 0, 1, 2, 3, 4, 5]);
const sources = {
  'resource-0001': { kind: 'web-file', name: 'genome.gb', type: 'text/plain', size: genome.byteLength, lastModified: 7 },
  'resource-0002': { name: 'reads.tsv.gz', size: gzipInput.byteLength, encoding: 'base64', data: 'unused', checksum: 'abc' }
};
const containerResources = sessionContainerResources(sources);
assert.deepEqual(containerResources, {
  'resource-0001': { kind: 'web-file', name: 'genome.gb', type: 'text/plain', size: genome.byteLength, lastModified: 7 },
  'resource-0002': { kind: 'web-file', name: 'reads.tsv.gz', type: 'application/octet-stream', size: gzipInput.byteLength, lastModified: 0 }
});
const payload = { format: 'gbdraw-session', version: 51, resources: containerResources };
const read = [];
const saved = await writeSessionFile(payload, {
  ...NO_COPIES,
  readResource: async (resourceId) => {
    read.push(resourceId);
    return streamOf(resourceId === 'resource-0001' ? genome : gzipInput);
  }
});
assert.deepEqual(read, ['resource-0001', 'resource-0002']);
assert.equal(saved.blob.type, 'application/zip');
const container = await readSessionContainer(saved.blob);
assert.deepEqual(container.document, payload);
assert.deepEqual(Object.keys(container.resources), ['resource-0001', 'resource-0002']);
assert.deepEqual([container.resources['resource-0001'].method, container.resources['resource-0002'].method], [8, 0]);
const entryBytes = async (resourceId) => new Uint8Array(await new Response(
  openSessionContainerEntry(saved.blob, container.resources[resourceId])).arrayBuffer());
assert.deepEqual(await entryBytes('resource-0001'), genome);
assert.deepEqual(await entryBytes('resource-0002'), gzipInput);

// Each written entry waits in a temporary copy when the browser gives one (J3); the
// file is read from the copies, which the caller gives back after the download.
const released = [];
const copyKeeper = {
  keep: async (data) => {
    const blob = await new Response(data).blob();
    const id = released.length + read.length;
    read.push(`copy-${id}`);
    return { copy: { blob: async () => blob, read: async () => new Uint8Array(await blob.arrayBuffer()), release: () => released.push(id) } };
  }
};
read.length = 0;
const kept = await writeSessionFile(payload, {
  ...copyKeeper, readResource: async (resourceId) => streamOf(resourceId === 'resource-0001' ? genome : gzipInput)
});
assert.equal(read.length, 3, 'the document and each input wait in a copy');
assert.deepEqual((await readSessionContainer(kept.blob)).document, payload);
assert.deepEqual(released, []);
kept.release();
assert.equal(released.length, 3);

// A copy store whose first write finds no storage: the entries wait in memory instead.
const memorySave = await writeSessionFile(payload, {
  keep: async () => ({ copy: null, reason: 'unavailable' }),
  readResource: async (resourceId) => streamOf(resourceId === 'resource-0001' ? genome : gzipInput)
});
assert.deepEqual((await readSessionContainer(memorySave.blob)).document, payload);

// A temporary write the browser refuses deletes what was written and refuses the
// Save with the storage text's facts: what it needs (inputs and settings and results),
// what is free, and the largest inputs; with no estimate, no free size.
const refusedSave = async (estimate) => {
  const releasedCopies = [];
  let entries = 0;
  const originalStorage = Object.getOwnPropertyDescriptor(globalThis.navigator, 'storage');
  Object.defineProperty(globalThis.navigator, 'storage', { configurable: true, value: estimate ? { estimate } : undefined });
  try {
    const error = await writeSessionFile({ ...payload, webFiles: { resourceOriginalNames: { 'resource-0001': 'My genome.gb' } } }, {
      readResource: async (resourceId) => streamOf(resourceId === 'resource-0001' ? genome : gzipInput),
      keep: async (data) => {
        entries += 1;
        await new Response(data).arrayBuffer();
        if (entries === 3) return { copy: null, reason: 'write-failed' };
        return { copy: { blob: async () => new Blob([]), read: async () => new Uint8Array(), release: () => releasedCopies.push(entries) } };
      }
    }).then(() => assert.fail('a refused temporary write must not save'), (caught) => caught);
    return { error, releasedCopies };
  } finally {
    if (originalStorage) Object.defineProperty(globalThis.navigator, 'storage', originalStorage);
    else delete globalThis.navigator.storage;
  }
};
const documentBytes = Buffer.byteLength(JSON.stringify({ ...payload, webFiles: { resourceOriginalNames: { 'resource-0001': 'My genome.gb' } } }));
const withEstimate = await refusedSave(async () => ({ quota: 5000, usage: 1000 }));
assert.deepEqual([withEstimate.error.code, withEstimate.error.stage, withEstimate.error.operation], ['SESSION_SIZE_LIMIT', 'write', 'session-save']);
assert.deepEqual(withEstimate.error.context, {
  storageNeedBytes: documentBytes + genome.byteLength + gzipInput.byteLength, storageFreeBytes: 4000,
  largestInputs: 'My genome.gb, reads.tsv.gz'
});
assert.equal(withEstimate.releasedCopies.length, 2, 'the copies written before the refusal are given back');
const withoutEstimate = await refusedSave(null);
assert.equal(withoutEstimate.error.context.storageFreeBytes, undefined);
assert.equal(withoutEstimate.error.context.storageNeedBytes, documentBytes + genome.byteLength + gzipInput.byteLength);

// A resource whose bytes are not its declared size is refused, not written short.
await assert.rejects(writeSessionFile(payload, { ...NO_COPIES, readResource: async () => streamOf(new Uint8Array(3)) }),
  (error) => error.code === 'RESOURCE_INVALID');

const semanticPayload = {
  resources: {},
  ordered: { zebra: 1, alpha: 2, nested: [true, null, { value: 'kept' }] },
  omitted: undefined,
  numbers: [Number.NaN, Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY, -0, 1.25],
  sparse: [, undefined, () => {}, Symbol('array-value')],
  escaped: 'quote:" backslash:\\ controls:\b\f\n\r\t\u0000',
  unicode: '日本語 😀 café',
  loneSurrogates: `left:\ud800 right:\udfff`
};
Object.defineProperty(semanticPayload, 'symbolValue', {
  enumerable: true,
  value: Symbol('object-value')
});
const semanticRoundTrip = await assertJsonParity(semanticPayload);
assert.deepEqual(Object.keys(semanticRoundTrip.ordered), ['zebra', 'alpha', 'nested']);
assert.equal(Object.hasOwn(semanticRoundTrip, 'omitted'), false);
assert.deepEqual(semanticRoundTrip.numbers, [null, null, null, 0, 1.25]);
assert.deepEqual(semanticRoundTrip.sparse, [null, null, null, null]);
assert.equal(Object.hasOwn(semanticRoundTrip, 'symbolValue'), false);

const sliceBoundary = Math.floor((256 * 1024) / 8);
const boundaryString = `${'x'.repeat(sliceBoundary - 1)}😀${'y'.repeat(sliceBoundary + 5)}`;
const longString = `${boundaryString}${'\u0001'.repeat(300_000)}${boundaryString}`;
assert.equal((await assertJsonParity({ resources: {}, longString })).longString, longString);

const immutablePayload = Object.freeze({
  resources: Object.freeze({}),
  top: Object.freeze({ value: 'unchanged' }),
  rows: Object.freeze([Object.freeze({ id: 1 })])
});
const immutableBefore = structuredClone(immutablePayload);
await assertJsonParity(immutablePayload);
assert.deepEqual(immutablePayload, immutableBefore);

let accessorReads = 0;
const accessorPayload = { resources: {}, rows: [{ id: 1 }, { id: 2 }] };
Object.defineProperty(accessorPayload, 'derived', {
  enumerable: true,
  get() {
    accessorReads += 1;
    return { value: 'read-once' };
  }
});
const accessorExpected = JSON.stringify(accessorPayload);
accessorReads = 0;
assert.equal(await documentText(accessorPayload), accessorExpected);
assert.equal(accessorReads, 1, 'bounded subtree probing must not invoke accessors twice');

// A Gallery Session's settings and results, as Save writes them.
const galleryContainer = new URL('../../gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw', import.meta.url);
const galleryDocument = existsSync(galleryContainer)
  ? (await readSessionContainer(new Blob([readFileSync(galleryContainer)]))).document
  : JSON.parse(readFileSync(new URL('../../gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json', import.meta.url), 'utf8'));
const currentSession = { ...galleryDocument, resources: {} };
assert.deepEqual(await assertJsonParity(currentSession), currentSession);

const cyclic = { resources: {} };
cyclic.self = cyclic;
await assert.rejects(writeDocument(cyclic), /cyclic value/);
await assert.rejects(writeDocument({ resources: {}, value: 1n }), /BigInt/);
await assert.rejects(writeDocument({ resources: {}, value: new Date() }), /unsupported object/);
await assert.rejects(writeDocument({ resources: {}, value: new Map() }), /unsupported object/);
await assert.rejects(
  writeDocument({ resources: {}, value: { toJSON: () => 'unsupported' } }),
  /unsupported object/
);

// Save refuses, before anything is written, settings and results
// (`gbdraw-session.json`) over the browser's cap (`MAX_EXPANDED_SESSION_BYTES`);
// the inputs have no cap. The cap is a test hook's small stand-in here; at or
// under it the bytes are those of JSON.stringify, as before.
assert.equal(MAX_EXPANDED_SESSION_BYTES, 512 * 1024 * 1024);
const sizedSession = { format: 'gbdraw-session', resources: {}, results: ['B'.repeat(700000), 'C'.repeat(400000)] };
const sizedBytes = Buffer.byteLength(JSON.stringify(sizedSession));
const withCap = async (cap, run) => {
  globalThis.__GBDRAW_TEST_HOOKS__ = { maxExpandedSessionBytes: cap };
  try {
    return await run();
  } finally {
    delete globalThis.__GBDRAW_TEST_HOOKS__;
  }
};
const atCap = await withCap(sizedBytes, () => documentText(sizedSession));
assert.equal(atCap, JSON.stringify(sizedSession), 'a Session exactly at the cap is written unchanged');
for (const keeper of [NO_COPIES, copyKeeper]) {
  const refusal = await withCap(sizedBytes - 1, () => writeDocument(sizedSession, keeper)).then(
    () => assert.fail('a Session over the cap must not be written'),
    (error) => error
  );
  assert.deepEqual([refusal.code, refusal.stage, refusal.operation], ['SESSION_SIZE_LIMIT', 'write', 'session-save']);
  assert.deepEqual(refusal.context, { sessionBytes: sizedBytes, sessionByteLimit: sizedBytes - 1 });
}
// A browser's body reader may reject with a new TypeError instead of the
// stream's error (Chromium does); the refusal reaches the caller as it was.
const NativeResponse = globalThis.Response;
globalThis.Response = class {
  constructor(body) { this.body = body; }
  async blob() {
    await new NativeResponse(this.body).blob().catch(() => {});
    throw new TypeError('Failed to fetch');
  }
};
try {
  const replaced = await withCap(sizedBytes - 1, () => writeDocument(sizedSession)).then(
    () => assert.fail('a Session over the cap must not be written'),
    (error) => error
  );
  assert.deepEqual([replaced.code, replaced.context.sessionBytes], ['SESSION_SIZE_LIMIT', sizedBytes]);
} finally {
  globalThis.Response = NativeResponse;
}
const blob = saved.blob;

const largeBlob = new Blob(
  [new Uint8Array((50 * 1024 * 1024) + (512 * 1024))],
  { type: 'application/gzip' }
);
// A Session over 50 MB asks first in the app's choice dialog, which
// names its size; a smaller one downloads without asking.
assert.equal(largeSessionDownloadMegabytes(largeBlob).toFixed(1), '50.5');
assert.equal(largeSessionDownloadMegabytes(new Blob([new Uint8Array(1024)])), null);
assert.throws(() => largeSessionDownloadMegabytes('not a blob'), TypeError);

// The question of every large download (Session, LOSAT tables, CLI helper
// files): one Download action; only Download downloads, Cancel and Escape
// (`cancel`) do not.
const asked = [];
const askAnswering = (answer) => async (request) => { asked.push(request); return answer; };
const question = { title: 'Download a large Session?', description: 'The Session file is 50.5 MB.' };
assert.equal(await askBeforeLargeDownload(askAnswering('download'), question), true);
assert.equal(await askBeforeLargeDownload(askAnswering('cancel'), question), false);
assert.deepEqual(asked[0], { ...question, actions: [{ id: 'download', label: 'Download', tone: 'primary' }] });

let downloadedBlob = null;
let downloadedName = null;
let clicked = false;
const originalCreateObjectUrl = URL.createObjectURL;
const originalRevokeObjectUrl = URL.revokeObjectURL;
URL.createObjectURL = (value) => {
  downloadedBlob = value;
  return 'blob:session-test';
};
URL.revokeObjectURL = () => {};
globalThis.document = {
  createElement: () => ({
    href: '',
    set download(value) {
      downloadedName = value;
    },
    click: () => {
      clicked = true;
    }
  })
};
try {
  downloadBlob(blob, 'exact-session.gbdraw');
  assert.equal(downloadedBlob, blob, 'the measured Session file is downloaded without copying');
  assert.equal(downloadedName, 'exact-session.gbdraw');
  assert.equal(clicked, true);
} finally {
  URL.createObjectURL = originalCreateObjectUrl;
  URL.revokeObjectURL = originalRevokeObjectUrl;
  delete globalThis.document;
}
