import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const {
  base64ToBytes,
  bytesToBase64,
  bytesToText,
  cloneFileBytesForTransfer,
  holdReadBytes,
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

// A single take hands out the bytes a read under the same hold cached, with no copy.
{
  const single = makeCountingFile('LOCUS single');
  const releaseReads = holdReadBytes();
  const cached = await readFileBytes(single);
  assert.equal(await takeFileBytesForTransfer(single), cached.buffer, 'the common take copies nothing');
  releaseReads();
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

// A gzip input expands to at most MAX_EXPANDED_INPUT_BYTES (Python's cap, from the
// generated contracts): past it the read rejects with INPUT_UNREADABLE, naming the
// file and the cap in MiB. A test hook stands in for the 512 MiB cap.
{
  const { MAX_EXPANDED_INPUT_BYTES } = await import('../../gbdraw/web/js/contracts/schema.generated.js');
  assert.equal(MAX_EXPANDED_INPUT_BYTES, 512 * 1024 * 1024);
  const cap = 1024 * 1024;
  const gzipOf = (name, size) => ({ name, async arrayBuffer() { return Uint8Array.from(gzipSync(Buffer.alloc(size, 0x41))).buffer; } });
  globalThis.__GBDRAW_TEST_HOOKS__ = { maxExpandedInputBytes: cap };
  try {
    assert.equal((await readFileBytes(gzipOf('at-cap.gb.gz', cap))).byteLength, cap);
    await assert.rejects(readFileBytes(gzipOf('big genome.gbff.gz', cap + 1)), (error) => {
      assert.deepEqual([error.code, error.stage, error.context],
        ['INPUT_UNREADABLE', 'read', { reason: 'GZIP_EXPANDED_SIZE', inputName: 'big genome.gbff.gz', limit: 1 }]);
      return true;
    });
    // A plain file larger than the cap reads whole.
    const plainBig = { name: 'plain.gb', async arrayBuffer() { return new Uint8Array(cap * 2).buffer; } };
    assert.equal((await readFileBytes(plainBig)).byteLength, cap * 2);
  } finally {
    delete globalThis.__GBDRAW_TEST_HOOKS__;
  }
}

// OV-2563: a File the browser refuses to read (Chromium rejects arrayBuffer()
// with NotReadableError for a file changed or removed on disk after it was
// picked) rejects at the read owner with INPUT_UNREADABLE, so every surface
// shows "An input could not be read" instead of the unknown failure. Another
// failure keeps its own error.
{
  const { normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');
  for (const name of ['NotReadableError', 'NotFoundError', 'SecurityError']) {
    const refused = { name: 'gone.gb', async arrayBuffer() { throw new DOMException('refused', name); } };
    const error = await readFileBytes(refused).then(() => null, (cause) => cause);
    assert.equal(error?.code, 'INPUT_UNREADABLE', name);
    assert.equal(error?.stage, 'read', name);
    const model = normalizeUserFacingError(error, { operation: 'input-records', stage: 'helper' });
    assert.equal(model?.code, 'INPUT_UNREADABLE', name);
    assert.equal(model?.summary, 'An input could not be read. Replace or reselect it.', name);
    const textError = await readFileText(refused).then(() => null, (cause) => cause);
    assert.equal(textError?.code, 'INPUT_UNREADABLE', `${name} through readFileText`);
  }
  const aborted = { name: 'other.gb', async arrayBuffer() { throw new DOMException('stopped', 'AbortError'); } };
  const other = await readFileBytes(aborted).then(() => null, (cause) => cause);
  assert.equal(other?.name, 'AbortError', 'a failure that is not a refused read keeps its own error');
}

// PD-OI-116, OV-2947, OV-2948: every input file the page reads through this
// module may be gzip (GenBank/DDBJ/EMBL, GFF3, FASTA, TSV tables, depth tables,
// comparison tables), so every input picker of the page offers gzip files. The
// Session picker reads Sessions, not input files, and the uploader template
// binds the accept list its user passes.
{
  const { readFile } = await import('node:fs/promises');
  const html = await readFile(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const pickers = [...html.matchAll(/\saccept="([^"]*)"/g)].map((match) => match[1])
    .filter((accept) => !accept.includes('.gbdraw') && !accept.startsWith('accept ||'));
  assert.ok(pickers.length >= 20, `found ${pickers.length} input pickers`);
  assert.deepEqual(pickers.filter((accept) => !accept.split(',').includes('.gz') || !accept.split(',').includes('application/gzip')), []);
}
