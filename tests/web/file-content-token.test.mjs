// The content token of an input and how long the page keeps the input's bytes
// (memory plan P2): one token scheme for every size, computed from the read the
// page already does, and the bytes dropped once no reader holds them. A File
// changed on disk after its bytes were dropped is refused at its next read (PD-OI-120).
import assert from 'node:assert/strict';
import test from 'node:test';
import { createHash, webcrypto } from 'node:crypto';
import { gzipSync } from 'node:zlib';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const {
  contentTokenOfBytes, holdReadBytes, readFileBytes, readFileContentFacts
} = await import('../../gbdraw/web/js/services/file-content-cache.js');
const { createDiagramResourceTransport } = await import('../../gbdraw/web/js/services/diagram-resource-staging.js');
const { setResourcePayloadOwner } = await import('../../gbdraw/web/js/services/resource-payload-owner.js');

const PIECE = 16 * 1024 * 1024;
const sha256 = (bytes) => createHash('sha256').update(bytes).digest();
// The agreed scheme: SHA-256 of the SHA-256 digests of the 16 MiB pieces, in order.
const expectedToken = (bytes) => {
  const digests = [];
  for (let start = 0; start === 0 || start < bytes.byteLength; start += PIECE) {
    digests.push(sha256(bytes.subarray(start, start + PIECE)));
  }
  return sha256(Buffer.concat(digests)).toString('hex');
};

/** A File whose every read is counted. */
const countingFile = (bytes, name = 'input.tsv') => ({
  name,
  reads: 0,
  async arrayBuffer() {
    this.reads += 1;
    return bytes.slice().buffer;
  }
});

test('the token is the digest of the 16 MiB piece digests, below, at, and above one piece', async () => {
  for (const size of [0, 5, PIECE, PIECE + 1, 2 * PIECE + 3]) {
    const bytes = new Uint8Array(size);
    for (let index = 0; index < size; index += 4093) bytes[index] = index % 251;
    const token = await contentTokenOfBytes(bytes);
    assert.match(token, /^[0-9a-f]{64}$/, `size ${size}`);
    assert.equal(token, expectedToken(Buffer.from(bytes.buffer)), `size ${size}`);
  }
});

test('the read that discovers a File gives its token and length, and the page then drops its bytes', async () => {
  const bytes = new TextEncoder().encode('LOCUS a\n//\n');
  const file = countingFile(bytes);
  await readFileBytes(file);
  const facts = await readFileContentFacts(file);
  assert.equal(file.reads, 1, 'the token comes from the read the page already did');
  assert.deepEqual(facts, { token: expectedToken(Buffer.from(bytes)), byteLength: bytes.byteLength });
  await readFileBytes(file);
  assert.equal(file.reads, 2, 'the bytes were dropped once the read settled');
  assert.equal((await readFileContentFacts(file)).token, facts.token);
  assert.equal(file.reads, 2, 'the token is never computed again');
});

test('reads that one hold covers share one read, and the bytes are dropped when the hold ends', async () => {
  const file = countingFile(new TextEncoder().encode('>x\nACGT\n'));
  const release = holdReadBytes();
  await readFileBytes(file);
  await readFileBytes(file);
  assert.equal(file.reads, 1);
  release();
  release();
  await readFileBytes(file);
  assert.equal(file.reads, 2);
});

test('a gzip input keeps its expanded bytes, and its token hashes them', async () => {
  const plain = new TextEncoder().encode('chr\t1\t5\n'.repeat(64));
  const file = countingFile(new Uint8Array(gzipSync(plain)), 'depth.tsv.gz');
  assert.deepEqual(await readFileBytes(file), plain);
  await readFileBytes(file);
  assert.equal(file.reads, 1, 'expanding again would be slower than keeping');
  assert.deepEqual(await readFileContentFacts(file), {
    token: expectedToken(Buffer.from(plain)), byteLength: plain.byteLength
  });
});

const fileRequest = { records: [{ source: { kind: 'genbank', resourceId: 'source' } }], diagramOptions: {} };
const fileResources = (file, size) => ({
  source: setResourcePayloadOwner({ kind: 'genbank', name: file.name, type: 'text/plain', size, lastModified: 1,
    encoding: 'base64', data: null }, file)
});

test('a File the browser refuses after its bytes were dropped fails its next read with INPUT_UNREADABLE', async () => {
  const bytes = new TextEncoder().encode('LOCUS c\n//\n');
  const file = countingFile(bytes);
  await readFileBytes(file);
  // The file changed on disk: Chrome refuses every later read of this File.
  file.arrayBuffer = async () => { throw new DOMException('changed', 'NotReadableError'); };
  await assert.rejects(readFileBytes(file), (error) => error.code === 'INPUT_UNREADABLE' && error.stage === 'read');
  assert.equal((await readFileContentFacts(file)).byteLength, bytes.byteLength, 'the facts from the first read stay');
});

test('Generate after the bytes were dropped reads the File again and sends the same bytes', async () => {
  const bytes = new TextEncoder().encode('LOCUS b\nORIGIN\n//\n');
  const file = countingFile(bytes);
  await readFileBytes(file);
  await readFileContentFacts(file);
  const prepared = await createDiagramResourceTransport().prepare({
    request: fileRequest, resources: fileResources(file, bytes.byteLength)
  });
  assert.equal(file.reads, 2);
  assert.deepEqual(new Uint8Array(prepared.stagedResources[0].bytes), bytes);
  assert.equal(prepared.resourceManifest[0].cacheToken, expectedToken(Buffer.from(bytes)));
});

test('a File whose token the Worker holds is read nothing again at Generate', async () => {
  const bytes = new TextEncoder().encode('LOCUS c\nORIGIN\n//\n');
  const transport = createDiagramResourceTransport();
  const earlier = countingFile(bytes);
  const first = await transport.prepare({ request: fileRequest, resources: fileResources(earlier, bytes.byteLength) });
  assert.equal(first.commit(), true);
  const chosen = countingFile(bytes);
  await readFileBytes(chosen);
  await readFileContentFacts(chosen);
  const second = await transport.prepare({ request: fileRequest, resources: fileResources(chosen, bytes.byteLength) });
  assert.equal(second.stagedResources.length, 0);
  assert.equal(chosen.reads, 1, 'only the read that discovered it');
});

test('a staging shares one read between a token it computes and the transfer', async () => {
  const bytes = new TextEncoder().encode('chr\t1\t7\n');
  const file = countingFile(bytes);
  const prepared = await createDiagramResourceTransport().prepare({
    request: fileRequest, resources: fileResources(file, bytes.byteLength)
  });
  assert.equal(file.reads, 1);
  assert.deepEqual(new Uint8Array(prepared.stagedResources[0].bytes), bytes);
});
