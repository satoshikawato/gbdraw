import assert from 'node:assert/strict';
import { readFile, writeFile, mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'utils', 'zip.js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-zip-'));
const modulePath = join(tempDir, 'zip.mjs');
await writeFile(modulePath, await readFile(sourcePath, 'utf8'), 'utf8');

const { createZipBlob } = await import(pathToFileURL(modulePath));

const readUint16 = (bytes, offset) => bytes[offset] | (bytes[offset + 1] << 8);
const readUint32 = (bytes, offset) =>
  (bytes[offset] |
    (bytes[offset + 1] << 8) |
    (bytes[offset + 2] << 16) |
    (bytes[offset + 3] << 24)) >>> 0;

const parseLocalFiles = (bytes) => {
  const files = [];
  const decoder = new TextDecoder();
  let offset = 0;
  while (readUint32(bytes, offset) === 0x04034B50) {
    const compressedSize = readUint32(bytes, offset + 18);
    const nameLength = readUint16(bytes, offset + 26);
    const extraLength = readUint16(bytes, offset + 28);
    const nameStart = offset + 30;
    const dataStart = nameStart + nameLength + extraLength;
    const name = decoder.decode(bytes.slice(nameStart, nameStart + nameLength));
    const data = bytes.slice(dataStart, dataStart + compressedSize);
    files.push({ name, text: decoder.decode(data), bytes: data });
    offset = dataStart + compressedSize;
  }
  return files;
};

const blob = createZipBlob([
  { name: 'alpha.tsv', data: 'a\tb\n' },
  { name: '../beta.tsv', data: 'b\tc\n' },
  { name: 'alpha.tsv', data: 'second\n' }
]);

assert.equal(blob.type, 'application/zip');
const bytes = new Uint8Array(await blob.arrayBuffer());
assert.equal(readUint32(bytes, 0), 0x04034B50);
assert.equal(readUint32(bytes, bytes.length - 22), 0x06054B50);
assert.equal(readUint16(bytes, bytes.length - 14), 3);

const files = parseLocalFiles(bytes);
assert.deepEqual(files.map((file) => file.name), ['alpha.tsv', 'beta.tsv', 'alpha-2.tsv']);
assert.deepEqual(files.map((file) => file.text), ['a\tb\n', 'b\tc\n', 'second\n']);

console.log('zip utility tests passed');
