import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';

const {
  MAX_EXPANDED_SESSION_BYTES,
  compressSessionData,
  largeSessionDownloadMegabytes
} = await import('../../gbdraw/web/js/services/session-file.js');
const { askBeforeLargeDownload, downloadBlob } = await import('../../gbdraw/web/js/services/text-download.js');

const compressedJsonText = async (value) => gunzipSync(
  Buffer.from(await (await compressSessionData(value)).arrayBuffer())
).toString('utf8');

const assertJsonParity = async (value) => {
  const actual = await compressedJsonText(value);
  assert.equal(actual, JSON.stringify(value));
  return JSON.parse(actual);
};

const payload = {
  format: 'gbdraw-session',
  version: 40,
  resources: {
    input: { data: 'A'.repeat(4096) }
  }
};
const blob = await compressSessionData(payload);
assert.equal(blob.type, 'application/gzip');
assert.equal(blob.size, blob.arrayBuffer ? (await blob.arrayBuffer()).byteLength : blob.size);
assert.deepEqual(
  JSON.parse(gunzipSync(Buffer.from(await blob.arrayBuffer())).toString('utf8')),
  payload
);

const semanticPayload = {
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
assert.equal((await assertJsonParity({ longString })).longString, longString);

const immutablePayload = Object.freeze({
  top: Object.freeze({ value: 'unchanged' }),
  rows: Object.freeze([Object.freeze({ id: 1 })])
});
const immutableBefore = structuredClone(immutablePayload);
await assertJsonParity(immutablePayload);
assert.deepEqual(immutablePayload, immutableBefore);

let accessorReads = 0;
const accessorPayload = { rows: [{ id: 1 }, { id: 2 }] };
Object.defineProperty(accessorPayload, 'derived', {
  enumerable: true,
  get() {
    accessorReads += 1;
    return { value: 'read-once' };
  }
});
const accessorExpected = JSON.stringify(accessorPayload);
accessorReads = 0;
assert.equal(await compressedJsonText(accessorPayload), accessorExpected);
assert.equal(accessorReads, 1, 'bounded subtree probing must not invoke accessors twice');

const currentSession = JSON.parse(readFileSync(new URL(
  '../../gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json',
  import.meta.url
), 'utf8'));
assert.deepEqual(await assertJsonParity(currentSession), currentSession);

const cyclic = {};
cyclic.self = cyclic;
await assert.rejects(compressSessionData(cyclic), /cyclic value/);
await assert.rejects(compressSessionData({ value: 1n }), /BigInt/);
await assert.rejects(compressSessionData({ value: new Date() }), /unsupported object/);
await assert.rejects(compressSessionData({ value: new Map() }), /unsupported object/);
await assert.rejects(
  compressSessionData({ value: { toJSON: () => 'unsupported' } }),
  /unsupported object/
);

// OV-2240: Save refuses, before anything is written, a Session whose expanded
// JSON is over the cap that Load enforces (`MAX_EXPANDED_SESSION_BYTES`). The
// cap is a test hook's small stand-in here; at or under it the bytes are those
// of JSON.stringify, as before.
assert.equal(MAX_EXPANDED_SESSION_BYTES, 512 * 1024 * 1024);
const sizedSession = {
  format: 'gbdraw-session',
  resources: {
    'resource-0001': { name: 'resource-0001-small.gb', size: 10, data: 'A'.repeat(10) },
    'resource-0002': { name: 'resource-0002-big_run_1_.tsv', size: 700000, data: 'B'.repeat(700000) },
    'resource-0003': { name: 'resource-0003-mid.fa', size: 400000, data: 'C'.repeat(400000) },
    'resource-0004': { name: 'resource-0004-tiny.tsv', size: 1, data: 'D' }
  },
  // The names the files were uploaded under; a resource without one shows its stored name.
  webFiles: { resourceOriginalNames: { 'resource-0002': 'big run (1).tsv', 'resource-0003': 'mid.fa', 'resource-0004': 'tiny.tsv' } }
};
const sizedBytes = Buffer.byteLength(JSON.stringify(sizedSession));
const withCap = async (cap, run, hook = 'maxExpandedSessionBytes') => {
  globalThis.__GBDRAW_TEST_HOOKS__ = { [hook]: cap };
  try {
    return await run();
  } finally {
    delete globalThis.__GBDRAW_TEST_HOOKS__;
  }
};
const atCap = await withCap(sizedBytes, () => compressSessionData(sizedSession));
assert.equal(
  gunzipSync(Buffer.from(await atCap.arrayBuffer())).toString('utf8'),
  JSON.stringify(sizedSession),
  'a Session exactly at the cap is written unchanged'
);
const refusal = await withCap(sizedBytes - 1, () => compressSessionData(sizedSession)).then(
  () => assert.fail('a Session over the cap must not be written'),
  (error) => error
);
assert.deepEqual(
  [refusal.code, refusal.stage, refusal.operation],
  ['SESSION_SIZE_LIMIT', 'write', 'session-save']
);
assert.deepEqual(refusal.context, {
  sessionBytes: sizedBytes,
  sessionByteLimit: sizedBytes - 1,
  largestInputs: 'big run (1).tsv, mid.fa, resource-0001-small.gb'
});
// OV-2245: Save also refuses a gzip file over the cap that Load reads
// (`MAX_SESSION_FILE_BYTES`), naming the file's size; a file exactly at the cap
// is written as before. The test hook stands in for 200 MiB.
const fileBytes = atCap.size;
const withFileCap = (cap, run) => withCap(cap, run, 'maxSessionFileBytes');
const atFileCap = await withFileCap(fileBytes, () => compressSessionData(sizedSession));
assert.equal(
  gunzipSync(Buffer.from(await atFileCap.arrayBuffer())).toString('utf8'),
  JSON.stringify(sizedSession),
  'a Session file exactly at the cap is written unchanged'
);
const fileRefusal = await withFileCap(fileBytes - 1, () => compressSessionData(sizedSession)).then(
  () => assert.fail('a Session file over the cap must not be written'),
  (error) => error
);
assert.deepEqual(
  [fileRefusal.code, fileRefusal.stage, fileRefusal.operation],
  ['SESSION_SIZE_LIMIT', 'write', 'session-save']
);
assert.deepEqual(fileRefusal.context, {
  sessionBytes: fileBytes,
  sessionByteLimit: fileBytes - 1,
  sessionSizeKind: 'file',
  largestInputs: 'big run (1).tsv, mid.fa, resource-0001-small.gb'
});
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
  const replaced = await withCap(sizedBytes - 1, () => compressSessionData(sizedSession)).then(
    () => assert.fail('a Session over the cap must not be written'),
    (error) => error
  );
  assert.deepEqual([replaced.code, replaced.context.sessionBytes], ['SESSION_SIZE_LIMIT', sizedBytes]);
} finally {
  globalThis.Response = NativeResponse;
}

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
  downloadBlob(blob, 'exact-session.json.gz');
  assert.equal(downloadedBlob, blob, 'the measured gzip Blob is downloaded without copying');
  assert.equal(downloadedName, 'exact-session.json.gz');
  assert.equal(clicked, true);
} finally {
  URL.createObjectURL = originalCreateObjectUrl;
  URL.revokeObjectURL = originalRevokeObjectUrl;
  delete globalThis.document;
}
