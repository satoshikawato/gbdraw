import assert from 'node:assert/strict';
import { gunzipSync } from 'node:zlib';

const {
  compressSessionData,
  confirmLargeSessionBlob,
  downloadSessionBlob
} = await import('../../gbdraw/web/js/services/session-file.js');

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
  downloadSessionBlob(blob, 'exact-session.json.gz');
  assert.equal(downloadedBlob, blob, 'the measured gzip Blob is downloaded without copying');
  assert.equal(downloadedName, 'exact-session.json.gz');
  assert.equal(clicked, true);
} finally {
  URL.createObjectURL = originalCreateObjectUrl;
  URL.revokeObjectURL = originalRevokeObjectUrl;
  delete globalThis.document;
}
