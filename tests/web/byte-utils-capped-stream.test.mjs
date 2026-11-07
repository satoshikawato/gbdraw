// The one capped reader of an expanded stream (services/byte-utils.js): a gzip
// input file (`expandGzipPieces`) and a Session (`readSessionText`) stop reading
// once the bytes pass their cap.
import assert from 'node:assert/strict';
import { gzipSync } from 'node:zlib';
import { test } from 'node:test';

const { expandGzipPieces, expandsAsGzip, joinBytes, readCappedStream } = await import('../../gbdraw/web/js/services/byte-utils.js');

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

// PD-OI-116: a gzip input expands up to the cap and is refused past it. Its
// expanded bytes arrive in pieces of the given size, each in its own buffer.
test('expandGzipPieces expands in pieces up to the cap and refuses past it', async () => {
  const content = new Uint8Array(4096).map((_, index) => index % 251);
  const gzip = Uint8Array.from(gzipSync(content));
  const overCap = () => new Error('gzip over the cap');
  const pieces = [];
  await expandGzipPieces(gzip, 4096, overCap, 1000, async (piece) => { pieces.push(piece); });
  assert.deepEqual(pieces.map((piece) => piece.byteLength), [1000, 1000, 1000, 1000, 96]);
  assert.ok(pieces.every((piece) => piece.byteOffset === 0 && piece.buffer.byteLength === piece.byteLength));
  assert.deepEqual(joinBytes(pieces), content);
  await assert.rejects(expandGzipPieces(gzip, 4095, overCap, 1000, async () => {}), /gzip over the cap/);
  const empty = [];
  await expandGzipPieces(Uint8Array.from(gzipSync(new Uint8Array(0))), 10, overCap, 1000, async (piece) => { empty.push(piece); });
  assert.deepEqual(empty.map((piece) => piece.byteLength), [0], 'an empty file is one empty piece');
  // Bytes that are not gzip are not expanded, whatever their size.
  assert.equal(expandsAsGzip(gzip), true);
  assert.equal(expandsAsGzip(new Uint8Array(8192)), false);
});
