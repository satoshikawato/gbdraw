// The Session file container (`<title>.gbdraw`, services/session-container.js):
// the Web reader and writer of the ZIP that Python owns (gbdraw/session_container.py).
// The shared vectors (tests/fixtures/session_container_vectors.json) pin the
// two to each other; a written container reads back entry by entry; Python's
// `zipfile` reads what the Web writes and the Web reads what Python writes
// (also with data descriptors and ZIP64); damaged files and a ZIP without
// gbdraw-session.json are refused with the receipt's texts; an entry is found
// by resource id only.
import assert from 'node:assert/strict';
import test from 'node:test';
import { execFileSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { gzipSync } from 'node:zlib';

const { isSessionContainerFile, openSessionContainerEntry, readSessionContainer, writeSessionContainer } = await import('../../gbdraw/web/js/services/session-container.js');
const { SESSION_DOCUMENT_ENTRY } = await import('../../gbdraw/web/js/contracts/schema.generated.js');
const { crc32Update, zipCentralRecord, zipEndRecords, zipLocalHeader } = await import('../../gbdraw/web/js/utils/zip.js');
const { normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');

const VECTORS = JSON.parse(readFileSync(new URL('../fixtures/session_container_vectors.json', import.meta.url), 'utf8'));
const encoder = new TextEncoder();
const DAMAGED = VECTORS.messages.damaged.web;
const NO_DOCUMENT = VECTORS.messages['no-document'].web;
const INVALID_JSON = VECTORS.messages['invalid-json'].web;

/** @param {Uint8Array} bytes */
const streamOf = (bytes) => new Blob([bytes]).stream();
/** @param {ReadableStream<Uint8Array>} stream */
const bytesOfStream = async (stream) => new Uint8Array(await new Response(stream).arrayBuffer());
/** @param {Uint8Array} bytes */
const sha256 = (bytes) => createHash('sha256').update(bytes).digest('hex');
/** @param {Blob} blob */
const bytesOfBlob = async (blob) => new Uint8Array(await blob.arrayBuffer());
/** @param {string} text */
const fromBase64 = (text) => new Uint8Array(Buffer.from(text, 'base64'));

const plain = encoder.encode('LOCUS       test 12 bp DNA\nORIGIN\n        1 acgtacgtacgt\n//\n'.repeat(400));
const gzipped = new Uint8Array(gzipSync(encoder.encode('>seq\nACGT\n'.repeat(1000))));
const descriptor = (name, size) => ({ kind: 'web-file', name, type: 'application/octet-stream', size, lastModified: 0 });
const sampleDocument = {
  format: 'gbdraw-session', version: 51,
  resources: { 'resource-1': descriptor('a.gb', plain.length), 'resource-2': descriptor('b.fa.gz', gzipped.length), 'resource-3': descriptor('c', 0) }
};
const documentText = JSON.stringify(sampleDocument);
const sampleBytes = { 'resource-1': plain, 'resource-2': gzipped, 'resource-3': new Uint8Array(0) };

/** Writes `document` with the bytes of `bytes` by id. */
const write = (document, bytes, options) => writeSessionContainer({
  document, documentBytes: () => streamOf(encoder.encode(JSON.stringify(document))), readResource: (id) => streamOf(bytes[id])
}, options);
const writeSample = (options) => write(sampleDocument, sampleBytes, options);

/** The bytes of a read container's resource. */
const resourceBytes = (blob, container, id) => bytesOfStream(openSessionContainerEntry(blob, container.resources[id]));

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
/** A document entry that declares resource `a` of `size` bytes. */
const documentDeclaring = (size, id = 'a') => ({
  name: SESSION_DOCUMENT_ENTRY, data: encoder.encode(JSON.stringify({ resources: { [id]: descriptor('a', size) } }))
});

/** What a reader reaches reading a vector's file: the document and each declared resource in full, or the refusal code. */
const readVector = async (blob) => {
  try {
    const container = await readSessionContainer(blob);
    const document = await bytesOfStream(openSessionContainerEntry(blob, container.documentEntry));
    const resources = {};
    for (const id of Object.keys(container.resources)) {
      const bytes = await resourceBytes(blob, container, id);
      resources[id] = { size: bytes.length, sha256: sha256(bytes) };
    }
    return { ok: true, document: { size: document.length, sha256: sha256(document) }, resources };
  } catch (error) {
    const { context, summary } = normalizeUserFacingError(error, { operation: 'session-load' });
    const code = { CONTAINER_DAMAGED: 'damaged', CONTAINER_NO_DOCUMENT: 'no-document', JSON_FORMAT: 'invalid-json' }[context.reason];
    assert.equal(summary, VECTORS.messages[code]?.web, `the text of ${code}`);
    return { ok: false, error: code };
  }
};

for (const vector of VECTORS.read) {
  test(`shared vector read ${vector.id}: ${vector.description}`, async () => {
    assert.deepEqual(await readVector(new Blob([fromBase64(vector.file)])), vector.expect);
  });
}

for (const vector of VECTORS.write) {
  test(`shared vector write ${vector.id}: ${vector.description}`, async () => {
    const bytes = Object.fromEntries(Object.entries(vector.resources).map(([id, text]) => [id, fromBase64(text)]));
    const blob = await write(vector.document, bytes);
    const file = await bytesOfBlob(blob);
    const view = new DataView(file.buffer);
    const entries = [];
    for (let at = 0; view.getUint32(at, true) === 0x04034B50;) {
      const nameLength = view.getUint16(at + 26, true);
      entries.push([new TextDecoder().decode(file.subarray(at + 30, at + 30 + nameLength)), view.getUint16(at + 8, true)]);
      at += 30 + nameLength + view.getUint16(at + 28, true) + view.getUint32(at + 18, true);
    }
    assert.deepEqual(entries, vector.expect.entries);
    const container = await readSessionContainer(blob);
    assert.deepEqual(container.document, vector.document);
    for (const [id, expected] of Object.entries(bytes)) assert.deepEqual(await resourceBytes(blob, container, id), expected);
  });
}


test('a written container reads back: the document first, each resource by its id, gzip bytes stored and others deflated', async () => {
  const blob = await writeSample();
  const container = await readSessionContainer(blob);
  assert.deepEqual(container.document, sampleDocument);
  assert.equal(container.resources['resource-1'].method, 8);
  assert.ok(container.resources['resource-1'].compressedSize < plain.length);
  assert.equal(container.resources['resource-2'].method, 0);
  assert.equal(container.resources['resource-3'].method, 8);
  for (const [id, bytes] of Object.entries(sampleBytes)) assert.deepEqual(await resourceBytes(blob, container, id), bytes);
  assert.equal(container.inputBytes, plain.length + gzipped.length);
  const bytes = await bytesOfBlob(blob);
  assert.equal(new TextDecoder().decode(bytes.subarray(30, 30 + SESSION_DOCUMENT_ENTRY.length)), SESSION_DOCUMENT_ENTRY);
  // An entry is plain data a Worker can hand to the page.
  assert.deepEqual(structuredClone(container.resources['resource-1']), container.resources['resource-1']);
});

test('a container written with every size in the ZIP64 fields reads back the same', async () => {
  const blob = await writeSample({ zip64: true });
  const bytes = await bytesOfBlob(blob);
  assert.equal(new DataView(bytes.buffer).getUint32(18, true), 0xFFFFFFFF, 'the local header points to the ZIP64 extra field');
  const container = await readSessionContainer(blob);
  for (const [id, expected] of Object.entries(sampleBytes)) assert.deepEqual(await resourceBytes(blob, container, id), expected);
});

test('each entry goes through the caller\'s keeper, and bytes of another size than declared are refused', async () => {
  const kept = [];
  await writeSample({ keep: async (data) => { const blob = await new Response(data).blob(); kept.push(blob.size); return blob; } });
  assert.equal(kept.length, 4);
  await assert.rejects(write({ resources: { a: descriptor('a', 5) } }, { a: plain }), (error) => error.code === 'RESOURCE_INVALID');
  await assert.rejects(write({ resources: { a: descriptor('a', -1) } }, { a: plain }), (error) => error.code === 'RESOURCE_INVALID');
});

test('the Web reader admits the document by its declared size before it reads it', async () => {
  const blob = await writeSample();
  const sizes = [];
  const tooLarge = new Error('too large');
  await assert.rejects(readSessionContainer(blob, { admitDocument: (size) => { sizes.push(size); throw tooLarge; } }), (error) => error === tooLarge);
  assert.deepEqual(sizes, [encoder.encode(documentText).length]);
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
        archive.writestr('resources/resource-3', b'')
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
    const blob = new Blob([bytes]);
    const container = await readSessionContainer(blob);
    assert.deepEqual(container.document, sampleDocument, tag);
    for (const [id, expected] of Object.entries(sampleBytes)) assert.deepEqual(await resourceBytes(blob, container, id), expected, tag);
  }
});

test('a file without a readable directory is refused as damaged: cut off, empty, not a ZIP, or a directory that does not add up', async () => {
  const bytes = await bytesOfBlob(await writeSample());
  await refused(readSessionContainer(new Blob([bytes.subarray(0, bytes.length - 10)])), DAMAGED);
  await refused(readSessionContainer(new Blob([bytes.subarray(0, bytes.length / 2)])), DAMAGED);
  await refused(readSessionContainer(new Blob([])), DAMAGED);
  await refused(readSessionContainer(new Blob([documentText])), DAMAGED);
  const end = bytes.length - 22;
  const wrongCount = bytes.slice();
  new DataView(wrongCount.buffer).setUint16(end + 10, 5, true);
  new DataView(wrongCount.buffer).setUint16(end + 8, 5, true);
  await refused(readSessionContainer(new Blob([wrongCount])), DAMAGED);
  const pastEnd = bytes.slice();
  new DataView(pastEnd.buffer).setUint32(end + 16, bytes.length, true);
  await refused(readSessionContainer(new Blob([pastEnd])), DAMAGED);
  // A directory offset below the real one reads as bytes before the ZIP (as in Python's zipfile): the
  // directory is found, and the document, looked for that much further on, is refused.
  const early = bytes.slice();
  new DataView(early.buffer).setUint32(end + 16, new DataView(bytes.buffer).getUint32(end + 16, true) - 7, true);
  await refused(readSessionContainer(new Blob([early])), DAMAGED);
});

test('a resource whose bytes, size or CRC-32 do not match is refused as damaged when it is read', async () => {
  const bytes = await bytesOfBlob(await writeSample());
  const storedStart = (() => { for (let i = 0; i < bytes.length - gzipped.length; i += 1) if (gzipped.every((b, k) => bytes[i + k] === b)) return i; return -1; })();
  assert.ok(storedStart > 0);
  const flipped = bytes.slice();
  flipped[storedStart + 10] ^= 0xFF;
  const flippedBlob = new Blob([flipped]);
  await refused(resourceBytes(flippedBlob, await readSessionContainer(flippedBlob), 'resource-2'), DAMAGED);
  const readA = async (blob) => resourceBytes(blob, await readSessionContainer(blob), 'a');
  // Deflated bytes that do not inflate.
  await refused(readA(craft([documentDeclaring(plain.length), { name: 'resources/a', data: plain, method: 8, compressed: new Uint8Array([0xFF, 0xFF, 0xFF, 0xFF]) }])), DAMAGED);
  // A size one byte larger, then one byte smaller, than the content, declared alike in the directory and the document.
  for (const size of [plain.length + 1, plain.length - 1]) {
    await refused(readA(craft([documentDeclaring(size), { name: 'resources/a', data: plain, size }])), DAMAGED);
  }
  // A CRC-32 in the directory that the content does not have.
  await refused(readA(craft([documentDeclaring(plain.length), { name: 'resources/a', data: plain, crc: 12345 }])), DAMAGED);
});

test('a declared resource without an entry, or whose entry has another size, is refused as damaged', async () => {
  await refused(readSessionContainer(craft([documentDeclaring(plain.length)])), DAMAGED);
  await refused(readSessionContainer(craft([documentDeclaring(plain.length + 1), { name: 'resources/a', data: plain }])), DAMAGED);
  await refused(readSessionContainer(craft([documentDeclaring(SESSION_DOCUMENT_ENTRY.length, '../gbdraw-session.json')])), DAMAGED);
});

test('two entries of one name, an encrypted entry and an unsupported compression are refused as damaged', async () => {
  const doc = documentDeclaring(plain.length);
  await refused(readSessionContainer(craft([doc, { name: 'resources/a', data: plain }, { name: 'resources/a', data: plain }])), DAMAGED);
  await refused(readSessionContainer(craft([doc, doc])), DAMAGED);
  const encrypted = (record) => { record[8] |= 0x01; };
  await refused(readSessionContainer(craft([doc, { name: 'resources/a', data: plain, central: encrypted }])), DAMAGED);
  const strong = (record) => { record[8] |= 0x40; };
  await refused(readSessionContainer(craft([doc, { name: 'resources/a', data: plain, central: strong }])), DAMAGED);
  await refused(readSessionContainer(craft([doc, { name: 'resources/a', data: plain, method: 12 }])), DAMAGED);
  // Any entry, also one the document does not declare.
  await refused(readSessionContainer(craft([doc, { name: 'resources/a', data: plain }, { name: 'notes.txt', data: plain, method: 14 }])), DAMAGED);
});

test('an entry whose local header names another entry, or whose data runs into the next entry, is refused as damaged', async () => {
  const doc = documentDeclaring(plain.length);
  const renamed = (header) => { header[30 + 'resources/'.length] = 'b'.charCodeAt(0); };
  const misnamed = craft([doc, { name: 'resources/a', data: plain, local: renamed }]);
  await refused(resourceBytes(misnamed, await readSessionContainer(misnamed), 'a'), DAMAGED);
  const longer = (record) => { new DataView(record.buffer).setUint32(20, plain.length + 50, true); };
  const overrun = craft([doc, { name: 'resources/a', data: plain, central: longer }, { name: 'resources/b', data: plain }]);
  await refused(resourceBytes(overrun, await readSessionContainer(overrun), 'a'), DAMAGED);
});

test('a ZIP without gbdraw-session.json is refused as no Session, and a document that is not UTF-8 JSON as not JSON', async () => {
  await refused(readSessionContainer(craft([{ name: 'resources/a', data: plain }])), NO_DOCUMENT);
  await refused(readSessionContainer(craft([{ name: 'gbdraw-session.json/', data: new Uint8Array(0) }])), NO_DOCUMENT);
  await refused(readSessionContainer(craft([{ name: SESSION_DOCUMENT_ENTRY, data: encoder.encode('{"resources": ') }])), INVALID_JSON);
  await refused(readSessionContainer(craft([{ name: SESSION_DOCUMENT_ENTRY, data: new Uint8Array([0x22, 0xFF, 0x22]) }])), INVALID_JSON);
});

// OV-2694 (D-30): Python refuses a document with one key twice as invalid JSON, so the Web does (until the shared
// vector case arrives, this local case pins it); one key in two objects, or a key-like text in a string, is no duplicate.
test('a document that holds one key twice is refused as not JSON, as Python refuses it', async () => {
  for (const text of ['{"resources": {}, "resources": {}}', '{"resources": {}, "a": {"b": 1, "\\u0062": 2}}']) {
    await refused(readSessionContainer(craft([{ name: SESSION_DOCUMENT_ENTRY, data: encoder.encode(text) }])), INVALID_JSON);
  }
  const kept = await readSessionContainer(craft([{ name: SESSION_DOCUMENT_ENTRY,
    data: encoder.encode('{"resources": {}, "a": [{"b": 1}, {"b": 2}], "c": "{\\"resources\\": 1}"}') }]));
  assert.deepEqual(kept.document, { resources: {}, a: [{ b: 1 }, { b: 2 }], c: '{"resources": 1}' });
});

test('content decides: a local header or an empty archive first is a container, any other file a JSON Session', async () => {
  const empty = craft([]);
  assert.equal(await isSessionContainerFile(empty), true);
  assert.equal(await isSessionContainerFile(craft([{ name: 'resources/a', data: plain }])), true);
  for (const bytes of ['{"format": "gbdraw-session"}', new Uint8Array([0x1f, 0x8b, 8, 0]), 'PK', '']) {
    assert.equal(await isSessionContainerFile(new Blob([bytes])), false);
  }
});

test('a resource is found by its id only: an entry name never becomes a path, and entries not declared are not read', async () => {
  const document = { resources: { '../x': descriptor('x', gzipped.length), '/x': descriptor('y', 1) } };
  const blob = craft([
    { name: SESSION_DOCUMENT_ENTRY, data: encoder.encode(JSON.stringify(document)) },
    { name: '../x', data: plain }, { name: 'resources/../x', data: gzipped }, { name: 'resources//x', data: encoder.encode('y') },
    { name: 'resources/x', data: plain, crc: 1 }
  ]);
  const container = await readSessionContainer(blob);
  assert.deepEqual(Object.keys(container.resources), ['../x', '/x']);
  assert.deepEqual(await resourceBytes(blob, container, '../x'), gzipped);
  assert.deepEqual(await resourceBytes(blob, container, '/x'), encoder.encode('y'));
  // Every entry but the document counts toward the storage a Load needs.
  assert.equal(container.inputBytes, 2 * plain.length + gzipped.length + 1);
});

test('a ZIP with bytes before it is read as Python\'s zipfile reads it', async () => {
  const bytes = await bytesOfBlob(await writeSample());
  const blob = new Blob([encoder.encode('#!prefix\n'), bytes]);
  assert.deepEqual(await resourceBytes(blob, await readSessionContainer(blob), 'resource-1'), plain);
});
