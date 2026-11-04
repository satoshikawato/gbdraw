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

// Each take is its caller's buffer to transfer: two helper requests that stage
// the same input at once (the region dialog derives on every input change,
// OV-2490) start their takes before either releases the cached bytes, and the
// second still receives bytes the first transfer leaves attached.
for (const owner of [rowA, makeCountingFile('LOCUS concurrent')]) {
  const [first, second] = await Promise.all([takeFileBytesForTransfer(owner), takeFileBytesForTransfer(owner)]);
  assert.notEqual(first, second, 'two takes never share one buffer');
  structuredClone(first, { transfer: [first] });
  assert.equal(first.byteLength, 0);
  assert.ok(second.byteLength > 0, 'the second take is not detached by the first transfer');
  structuredClone(second, { transfer: [second] });
}

// A single take hands out the cached buffer itself, with no copy.
{
  const single = makeCountingFile('LOCUS single');
  const cached = await readFileBytes(single);
  assert.equal(await takeFileBytesForTransfer(single), cached.buffer, 'the common take copies nothing');
}
// The retake is bounded: a FileLike that returns one stored buffer on every
// read gives the second take a copy instead of reading forever.
{
  const stored = encoder.encode('LOCUS stored').buffer;
  const sameBuffer = { name: 'same.txt', arrayBuffer: async () => stored };
  const [first, second] = await Promise.all([takeFileBytesForTransfer(sameBuffer), takeFileBytesForTransfer(sameBuffer)]);
  assert.equal(first, stored);
  assert.notEqual(second, stored, 'the second take is a copy');
  assert.equal(bytesToText(new Uint8Array(second)), 'LOCUS stored');
}

// A gzip input file (NCBI's .gbff.gz, .fna.gz, .gff.gz) reads as the plain
// file it holds, told by its first two bytes, not its name: every browser
// reader, the Session and the Worker see the plain bytes.
const { gzipSync } = await import('node:zlib');
const gzipFile = {
  name: 'input.gb.gz',
  async arrayBuffer() {
    return Uint8Array.from(gzipSync('LOCUS gzipped')).buffer;
  }
};
assert.equal(bytesToText(await readFileBytes(gzipFile)), 'LOCUS gzipped');
assert.equal(await readFileText(gzipFile), 'LOCUS gzipped');
const plainNamedGz = { name: 'plain.gb.gz', async arrayBuffer() { return encoder.encode('LOCUS plain').buffer; } };
assert.equal(await readFileText(plainNamedGz), 'LOCUS plain');
