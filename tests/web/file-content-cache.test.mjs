import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const {
  bytesToBase64,
  cloneFileBytesForTransfer,
  filePayloadIdentity,
  readFileBytes,
  readFileText
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
const [bytes, text, transfer, identity] = await Promise.all([
  readFileBytes(file),
  readFileText(file),
  cloneFileBytesForTransfer(file),
  filePayloadIdentity(file)
]);
assert.equal(file.reads, 1);
assert.equal(new TextDecoder().decode(bytes), 'cached bytes');
assert.equal(text, 'cached bytes');
assert.equal(new TextDecoder().decode(transfer), 'cached bytes');
assert.equal(identity.byteLength, 12);
assert.equal(identity.sha256, '6e10e7392858b23a91ec0fc838bdea89407aa958ff11f6d89dd6f593be5cc008');
assert.equal(bytesToBase64(bytes), 'Y2FjaGVkIGJ5dGVz');

const retry = makeCountingFile('retry', { failOnce: true });
await assert.rejects(readFileBytes(retry), /read failed/);
assert.equal(await readFileText(retry), 'retry');
assert.equal(retry.reads, 2);

const separate = makeCountingFile('cached bytes');
await readFileBytes(separate);
assert.equal(separate.reads, 1);
assert.equal(file.reads, 1);
