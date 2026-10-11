import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const {
  base64ToBytes,
  bytesToBase64,
  bytesToText,
  cloneFileBytesForTransfer,
  readFileBytes,
  readFileText,
  takeFileBytesForTransfer,
  textToBase64,
  textToBytes
} = await import('../../gbdraw/web/js/services/file-content-cache.js');

const encoder = new TextEncoder();
const makeCountingFile = (text, { failOnce = false } = {}) => {
  let reads = 0;
  return {
    name: 'input.txt',
    get reads() {
      return reads;
    },
    async arrayBuffer() {
      reads += 1;
      if (failOnce && reads === 1) throw new Error('read failed');
      return encoder.encode(text).buffer;
    }
  };
};

const file = makeCountingFile('cached bytes');
const [bytes, text, transfer] = await Promise.all([
  readFileBytes(file),
  readFileText(file),
  cloneFileBytesForTransfer(file)
]);
assert.equal(file.reads, 1);
assert.equal(new TextDecoder().decode(bytes), 'cached bytes');
assert.equal(text, 'cached bytes');
assert.equal(new TextDecoder().decode(transfer), 'cached bytes');
assert.equal(bytesToBase64(bytes), 'Y2FjaGVkIGJ5dGVz');
assert.equal(textToBase64('雪'), '6Zuq');
assert.equal(bytesToText(base64ToBytes('6Zuq')), '雪');
assert.deepEqual(textToBytes('AB'), new Uint8Array([65, 66]));

const retry = makeCountingFile('retry', { failOnce: true });
await assert.rejects(readFileBytes(retry), /read failed/);
assert.equal(await readFileText(retry), 'retry');
assert.equal(retry.reads, 2);

const separate = makeCountingFile('cached bytes');
await readFileBytes(separate);
assert.equal(separate.reads, 1);
assert.equal(file.reads, 1);

// A one-file two-record Linear Session fills two rows with File views of the
// same resource; the views share one backing. Transferring the bytes to the
// worker through one view detaches the shared buffer, and the other view reads
// the resource again instead of the detached array (OV-304).
const { adoptCurrentSessionResources, createSessionResourceFileView } = await import(
  '../../gbdraw/web/js/services/session-resource-backing.js'
);
const sharedTable = adoptCurrentSessionResources({
  'record-1-genbank': {
    kind: 'genbank', name: 'shared.gb', type: 'text/plain', size: 12, lastModified: 0,
    encoding: 'base64', data: textToBase64('LOCUS shared')
  }
});
const rowA = createSessionResourceFileView(sharedTable, 'record-1-genbank');
const rowB = createSessionResourceFileView(sharedTable, 'record-1-genbank');
assert.equal(bytesToText(await readFileBytes(rowA)), 'LOCUS shared');
assert.equal(bytesToText(await readFileBytes(rowB)), 'LOCUS shared');
const transferred = await takeFileBytesForTransfer(rowB);
structuredClone(transferred, { transfer: [transferred] });
assert.equal(transferred.byteLength, 0);
assert.equal(bytesToText(await readFileBytes(rowA)), 'LOCUS shared');
assert.equal(bytesToText(await readFileBytes(rowB)), 'LOCUS shared');
