// The one capped reader of an expanded stream (services/byte-utils.js): a gzip
// input file (`expandGzipBytes`) and a Session (`readSessionText`) stop reading
// once the bytes pass their cap.
import assert from 'node:assert/strict';
import { gzipSync } from 'node:zlib';
import { test } from 'node:test';

const { expandGzipBytes, readCappedStream } = await import('../../gbdraw/web/js/services/byte-utils.js');

/** @param {number[]} sizes */
const chunkStream = (sizes) => {
  const state = { cancelled: false, pulled: 0 };
  const stream = new ReadableStream({
    pull(controller) {
      if (state.pulled === sizes.length) return controller.close();
      controller.enqueue(new Uint8Array(sizes[state.pulled++]));
    },
    cancel() { state.cancelled = true; }
  });
  return { stream, state };
};

test('the capped reader passes every chunk up to the cap', async () => {
  const { stream } = chunkStream([3, 4, 3]);
  const seen = [];
  await readCappedStream(stream, 10, (chunk) => seen.push(chunk.byteLength), () => new Error('over'));
  assert.deepEqual(seen, [3, 4, 3]);
});

test('the capped reader stops past the cap with its error and cancels the stream', async () => {
  const { stream, state } = chunkStream([6, 6, 6]);
  const seen = [];
  await assert.rejects(
    readCappedStream(stream, 10, (chunk) => seen.push(chunk.byteLength), () => new Error('over the cap')),
    /over the cap/
  );
  assert.deepEqual(seen, [6], 'the chunk that passes the cap is not handed on');
  assert.equal(state.cancelled, true);
});

// PD-OI-116: a gzip input expands up to the cap and is refused past it.
test('expandGzipBytes expands up to the cap and refuses past it', async () => {
  const content = new Uint8Array(4096).fill(0x61);
  const gzip = Uint8Array.from(gzipSync(content));
  const overCap = () => new Error('gzip over the cap');
  assert.deepEqual(await expandGzipBytes(gzip, 4096, overCap), content);
  await assert.rejects(expandGzipBytes(gzip, 4095, overCap), /gzip over the cap/);
  // Bytes that are not gzip are returned as they are, whatever their size.
  const plain = new Uint8Array(8192);
  assert.equal(await expandGzipBytes(plain, 4095, overCap), plain);
});
