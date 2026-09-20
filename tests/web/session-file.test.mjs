import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';

const {
  compressSessionData,
  confirmLargeSessionBlob
} = await import('../../gbdraw/web/js/services/session-file.js');
const { downloadBlob } = await import('../../gbdraw/web/js/services/text-download.js');

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

const largeBlob = new Blob(
  [new Uint8Array((50 * 1024 * 1024) + (512 * 1024))],
  { type: 'application/gzip' }
);
let confirmationMessage = '';
assert.equal(
  confirmLargeSessionBlob(largeBlob, (message) => {
    confirmationMessage = message;
    return false;
  }),
  false
);
assert.equal(
  confirmationMessage,
  'Compressed session size is 50.5 MB. Continue?'
);
assert.equal(
  confirmLargeSessionBlob(largeBlob, () => true),
  true
);

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
