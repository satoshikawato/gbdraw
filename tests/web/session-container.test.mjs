// The Session file container (`<title>.gbdraw`, services/session-container.js):
// the Web reader and writer of the ZIP that Python owns (gbdraw/session_container.py).
// A written container reads back entry by entry; Python's `zipfile` reads what
// the Web writes and the Web reads what Python writes (also with data
// descriptors and ZIP64); damaged files and a ZIP without gbdraw-session.json
// are refused with the receipt's texts; an entry is found by resource id only.
import assert from 'node:assert/strict';
import test from 'node:test';
import { execFileSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { gzipSync } from 'node:zlib';

const { readSessionContainer, writeSessionContainer, SESSION_DOCUMENT_ENTRY } = await import('../../gbdraw/web/js/services/session-container.js');
const { crc32Update, zipCentralRecord, zipEndRecords, zipLocalHeader } = await import('../../gbdraw/web/js/utils/zip.js');
const { normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');

const encoder = new TextEncoder();
const DAMAGED = 'This file is not a gbdraw Session: it is incomplete or damaged.';
const NO_DOCUMENT = 'This file is not a gbdraw Session: it has no gbdraw-session.json.';

/** @param {Uint8Array} bytes */
const streamOf = (bytes) => () => new Blob([bytes]).stream();
/** @param {ReadableStream<Uint8Array>} stream */
const bytesOfStream = async (stream) => new Uint8Array(await new Response(stream).arrayBuffer());
/** @param {Uint8Array} bytes */
const sha256 = (bytes) => createHash('sha256').update(bytes).digest('hex');
/** @param {Blob} blob */
const bytesOfBlob = async (blob) => new Uint8Array(await blob.arrayBuffer());

const documentText = JSON.stringify({ format: 'gbdraw-session', version: 51, resources: {} });
const plain = encoder.encode('LOCUS       test 12 bp DNA\nORIGIN\n        1 acgtacgtacgt\n//\n'.repeat(400));
const gzipped = new Uint8Array(gzipSync(encoder.encode('>seq\nACGT\n'.repeat(1000))));

const writeSample = (options) => writeSessionContainer({
  document: streamOf(encoder.encode(documentText)),
  resources: [
    { id: 'resource-1', read: streamOf(plain), stored: false },
    { id: 'resource-2', read: streamOf(gzipped), stored: true },
    { id: 'resource-3', read: streamOf(new Uint8Array(0)), stored: false }
  ]
}, options);

/** @param {Promise<unknown>} promise @param {string} summary */
const refused = async (promise, summary) => {
  const error = await promise.then(() => assert.fail('expected a refusal'), (caught) => caught);
  assert.equal(normalizeUserFacingError(error, { operation: 'session-load' }).summary, summary);
};

/**
 * A ZIP built entry by entry, for the damaged files: each entry may change its
 * local header or central record after they are encoded.
 * @param {{ name: string, data: Uint8Array, method?: number, crc?: number, size?: number, local?: (header: Uint8Array) => void,
 *   central?: (record: Uint8Array) => void, compressed?: Uint8Array }[]} entries
 */
const craft = (entries) => {
  const parts = [];
  const directory = [];
  let offset = 0;
  for (const item of entries) {
    const data = item.compressed ?? item.data;
    const fields = {
      name: encoder.encode(item.name), method: item.method ?? 0, crc: item.crc ?? crc32Update(0, item.data),
      compressedSize: data.length, size: item.size ?? item.data.length, offset, time: 0, date: 33
    };
    const header = zipLocalHeader(fields);
    item.local?.(header);
    const record = zipCentralRecord(fields);
    item.central?.(record);
    parts.push(header, data);
    directory.push(record);
    offset += header.length + data.length;
  }
  const directorySize = directory.reduce((total, record) => total + record.length, 0);
  return new Blob([...parts, ...directory, zipEndRecords({ count: directory.length, directoryOffset: offset, directorySize })]);
};
const documentEntry = { name: SESSION_DOCUMENT_ENTRY, data: encoder.encode(documentText) };

test('a written container reads back: the document, each resource by its id, stored and deflated', async () => {
  const blob = await writeSample();
  const container = await readSessionContainer(blob);
  assert.equal(new TextDecoder().decode(await bytesOfStream(container.document.open())), documentText);
  const first = /** @type {NonNullable<ReturnType<typeof container.resource>>} */ (container.resource('resource-1'));
  assert.equal(first.stored, false);
  assert.ok(first.compressedSize < plain.length);
  assert.deepEqual(await bytesOfStream(first.open()), plain);
  const second = /** @type {NonNullable<ReturnType<typeof container.resource>>} */ (container.resource('resource-2'));
  assert.equal(second.stored, true);
  assert.equal(second.compressedSize, gzipped.length);
  assert.deepEqual(await bytesOfStream(second.open()), gzipped);
  assert.deepEqual(await bytesOfStream(container.resource('resource-3').open()), new Uint8Array(0));
  assert.equal(container.resource('resource-4'), null);
  assert.equal(container.resource(''), null);
  assert.equal(container.inputBytes, plain.length + gzipped.length);
  // The document is the first entry.
  const bytes = await bytesOfBlob(blob);
  assert.equal(new TextDecoder().decode(bytes.subarray(30, 30 + SESSION_DOCUMENT_ENTRY.length)), SESSION_DOCUMENT_ENTRY);
});

test('a container written with every size in the ZIP64 fields reads back the same', async () => {
  const blob = await writeSample({ zip64: true });
  const bytes = await bytesOfBlob(blob);
  assert.equal(new DataView(bytes.buffer).getUint32(18, true), 0xFFFFFFFF, 'the local header points to the ZIP64 extra field');
  const container = await readSessionContainer(blob);
  assert.deepEqual(await bytesOfStream(container.resource('resource-1').open()), plain);
  assert.deepEqual(await bytesOfStream(container.resource('resource-2').open()), gzipped);
});

test('an entry is kept by the caller\'s keeper, and two inputs of one id are refused before anything is written', async () => {
  const kept = [];
  await writeSample({ keep: async (data) => { const blob = await new Response(data).blob(); kept.push(blob.size); return blob; } });
  assert.equal(kept.length, 4);
  await assert.rejects(writeSessionContainer({
    document: streamOf(encoder.encode('{}')),
    resources: [{ id: 'a', read: streamOf(plain), stored: false }, { id: 'a', read: streamOf(plain), stored: false }]
  }), (error) => error.code === 'RESOURCE_INVALID');
});

test('Python\'s zipfile reads a container the Web wrote, and the Web reads one Python wrote, with data descriptors and ZIP64', async () => {
  const written = {
    plain: Buffer.from(await bytesOfBlob(await writeSample())).toString('base64'),
    zip64: Buffer.from(await bytesOfBlob(await writeSample({ zip64: true }))).toString('base64')
  };
  const report = JSON.parse(execFileSync('python', ['-c', `
import base64, hashlib, io, json, sys, zipfile
given = json.load(sys.stdin)
read = {}
for tag, text in given['written'].items():
    with zipfile.ZipFile(io.BytesIO(base64.b64decode(text))) as archive:
        assert archive.testzip() is None
        read[tag] = [[info.filename, info.compress_type, hashlib.sha256(archive.read(info)).hexdigest()] for info in archive.infolist()]
class Unseekable(io.RawIOBase):
    def __init__(self): self.parts = []
    def writable(self): return True
    def write(self, data): self.parts.append(bytes(data)); return len(data)
def build(target, zip64):
    with zipfile.ZipFile(target, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr('gbdraw-session.json', given['document'])
        with archive.open('resources/resource-1', 'w', force_zip64=zip64) as entry:
            entry.write(base64.b64decode(given['plain']))
        info = zipfile.ZipInfo('resources/resource-2')
        info.compress_type = zipfile.ZIP_STORED
        archive.writestr(info, base64.b64decode(given['gzipped']))
seekable = io.BytesIO(); build(seekable, False)
streamed = Unseekable(); build(streamed, True)
print(json.dumps({'read': read, 'seekable': base64.b64encode(seekable.getvalue()).decode(),
                  'streamed': base64.b64encode(b''.join(streamed.parts)).decode()}))
`], {
    input: JSON.stringify({
      written, document: documentText, plain: Buffer.from(plain).toString('base64'), gzipped: Buffer.from(gzipped).toString('base64')
    }),
    encoding: 'utf8', maxBuffer: 64 * 1024 * 1024
  }));
  const expected = [
    [SESSION_DOCUMENT_ENTRY, 8, sha256(encoder.encode(documentText))],
    ['resources/resource-1', 8, sha256(plain)], ['resources/resource-2', 0, sha256(gzipped)], ['resources/resource-3', 8, sha256(new Uint8Array(0))]
  ];
  assert.deepEqual(report.read.plain, expected);
  assert.deepEqual(report.read.zip64, expected);
  for (const tag of ['seekable', 'streamed']) {
    const bytes = new Uint8Array(Buffer.from(report[tag], 'base64'));
    if (tag === 'streamed') assert.equal(new DataView(bytes.buffer).getUint16(6, true) & 0x08, 0x08, 'Python wrote data descriptors');
    const container = await readSessionContainer(new Blob([bytes]));
    assert.equal(new TextDecoder().decode(await bytesOfStream(container.document.open())), documentText, tag);
    assert.deepEqual(await bytesOfStream(container.resource('resource-1').open()), plain, tag);
    assert.deepEqual(await bytesOfStream(container.resource('resource-2').open()), gzipped, tag);
    assert.equal(container.resource('resource-2').stored, true, tag);
  }
});

test('a file without a readable directory is refused as damaged: cut off, empty, not a ZIP, or a directory that does not add up', async () => {
  const bytes = await bytesOfBlob(await writeSample());
  await refused(readSessionContainer(new Blob([bytes.subarray(0, bytes.length - 10)])), DAMAGED);
  await refused(readSessionContainer(new Blob([bytes.subarray(0, bytes.length / 2)])), DAMAGED);
  await refused(readSessionContainer(new Blob([])), DAMAGED);
  await refused(readSessionContainer(new Blob([documentText])), DAMAGED);
  const wrongCount = bytes.slice();
  const end = wrongCount.length - 22;
  new DataView(wrongCount.buffer).setUint16(end + 10, 5, true);
  new DataView(wrongCount.buffer).setUint16(end + 8, 5, true);
  await refused(readSessionContainer(new Blob([wrongCount])), DAMAGED);
  const pastEnd = bytes.slice();
  new DataView(pastEnd.buffer).setUint32(end + 16, bytes.length, true);
  await refused(readSessionContainer(new Blob([pastEnd])), DAMAGED);
  // A directory offset below the real one reads as bytes before the ZIP (as in Python's zipfile): the
  // directory is found, and each entry, looked for that much further on, is refused when it is read.
  const early = bytes.slice();
  new DataView(early.buffer).setUint32(end + 16, new DataView(bytes.buffer).getUint32(end + 16, true) - 7, true);
  await refused(bytesOfStream((await readSessionContainer(new Blob([early]))).document.open()), DAMAGED);
});

test('an entry whose bytes, size or CRC-32 do not match is refused as damaged when it is read', async () => {
  const bytes = await bytesOfBlob(await writeSample());
  // A flipped byte in the stored entry: its CRC-32 differs.
  const storedStart = (() => { for (let i = 0; i < bytes.length - gzipped.length; i += 1) if (gzipped.every((b, k) => bytes[i + k] === b)) return i; return -1; })();
  assert.ok(storedStart > 0);
  const flipped = bytes.slice();
  flipped[storedStart + 10] ^= 0xFF;
  await refused(bytesOfStream((await readSessionContainer(new Blob([flipped]))).resource('resource-2').open()), DAMAGED);
  // Deflated bytes that do not inflate.
  const broken = craft([documentEntry, { name: 'resources/a', data: plain, method: 8, compressed: new Uint8Array([0xFF, 0xFF, 0xFF, 0xFF]) }]);
  await refused(bytesOfStream((await readSessionContainer(broken)).resource('a').open()), DAMAGED);
  // A declared size one byte larger, then one byte smaller, than the content.
  for (const size of [plain.length + 1, plain.length - 1]) {
    const wrongSize = craft([documentEntry, { name: 'resources/a', data: plain, size }]);
    const entry = (await readSessionContainer(wrongSize)).resource('a');
    await refused(bytesOfStream(entry.open()), DAMAGED);
  }
  // A CRC-32 in the directory that the content does not have.
  const wrongCrc = craft([documentEntry, { name: 'resources/a', data: plain, crc: 12345 }]);
  await refused(bytesOfStream((await readSessionContainer(wrongCrc)).resource('a').open()), DAMAGED);
});

test('two entries of one name, an encrypted entry and an unsupported compression are refused as damaged', async () => {
  await refused(readSessionContainer(craft([documentEntry, { name: 'resources/a', data: plain }, { name: 'resources/a', data: plain }])), DAMAGED);
  await refused(readSessionContainer(craft([documentEntry, documentEntry])), DAMAGED);
  const encrypted = (record) => { record[8] |= 0x01; };
  await refused(readSessionContainer(craft([documentEntry, { name: 'resources/a', data: plain, central: encrypted }])), DAMAGED);
  await refused(readSessionContainer(craft([documentEntry, { name: 'resources/a', data: plain, method: 12 }])), DAMAGED);
  await refused(readSessionContainer(craft([documentEntry, { name: 'resources/a', data: plain, method: 14 }])), DAMAGED);
});

test('an entry whose local header names another entry, or whose data runs into the next entry, is refused as damaged', async () => {
  const renamed = (header) => { header[30 + 'resources/'.length] = 'b'.charCodeAt(0); };
  const misnamed = await readSessionContainer(craft([documentEntry, { name: 'resources/a', data: plain, local: renamed }]));
  await refused(bytesOfStream(misnamed.resource('a').open()), DAMAGED);
  const longer = (record) => { new DataView(record.buffer).setUint32(20, plain.length + 50, true); };
  const overrun = await readSessionContainer(craft([documentEntry, { name: 'resources/a', data: plain, central: longer }, { name: 'resources/b', data: plain }]));
  await refused(bytesOfStream(overrun.resource('a').open()), DAMAGED);
});

test('a ZIP without gbdraw-session.json is refused as no Session', async () => {
  await refused(readSessionContainer(craft([{ name: 'resources/a', data: plain }])), NO_DOCUMENT);
  await refused(readSessionContainer(craft([{ name: 'gbdraw-session.json/', data: new Uint8Array(0) }])), NO_DOCUMENT);
});

test('an entry is found by its resource id only: a name never becomes a path', async () => {
  const container = await readSessionContainer(craft([
    documentEntry, { name: '../x', data: plain }, { name: 'resources/../x', data: gzipped }, { name: 'resources//x', data: encoder.encode('y') }
  ]));
  assert.equal(container.resource('x'), null);
  assert.deepEqual(await bytesOfStream(container.resource('../x').open()), gzipped);
  assert.deepEqual(await bytesOfStream(container.resource('/x').open()), encoder.encode('y'));
  assert.equal(container.resource('../../x'), null);
  // Every entry but the document counts toward the storage a Load needs.
  assert.equal(container.inputBytes, plain.length + gzipped.length + 1);
});

test('a ZIP with bytes before it is read as Python\'s zipfile reads it', async () => {
  const bytes = await bytesOfBlob(await writeSample());
  const container = await readSessionContainer(new Blob([encoder.encode('#!prefix\n'), bytes]));
  assert.deepEqual(await bytesOfStream(container.resource('resource-1').open()), plain);
});
