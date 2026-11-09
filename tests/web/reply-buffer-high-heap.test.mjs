// A Generate reply crosses from Python to the Worker whatever its address in the
// Pyodide heap (OV-2800). Pyodide 0.29 reads a Python buffer's address as a
// signed 32-bit number, so a chunk above 2 GiB failed with a RangeError when
// the Worker read it through a buffer view; it now crosses through a file.
// The real vendored Pyodide runs here (helpers/pyodide-high-heap.mjs).
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { HIGH_ADDRESS as HIGH, allocateAboveTwoGiB, loadVendoredPyodide } from './helpers/pyodide-high-heap.mjs';

globalThis.self = {};

const root = fileURLToPath(new URL('../../', import.meta.url));
const { createReplyBuffer, finishReplyBuffers } = await import('../../gbdraw/web/js/workers/diagram-generation-worker.js');

const pyodide = await loadVendoredPyodide();
pyodide.FS.writeFile('/home/pyodide/reply_buffers.py', readFileSync(`${root}gbdraw/web_support/reply_buffers.py`));

const release = allocateAboveTwoGiB(pyodide);

try {
  const write = pyodide.runPython(`
import hashlib
from reply_buffers import json_utf8_chunks, write_reply_buffer

def write(new_reply_buffer):
    value = {"items": [{"features": [{"product": "\\u00b5-protein \\u65e5\\u672c \\U0001F9EC", "n": i} for i in range(60000)]}]}
    digest = hashlib.sha256()
    addresses = []
    def chunks():
        for chunk in json_utf8_chunks(value):
            addresses.append(id(chunk))
            digest.update(chunk)
            yield chunk
    buffer = write_reply_buffer(new_reply_buffer, chunks())
    return {"buffer": buffer, "low": sum(address < ${HIGH} for address in addresses), "chunks": len(addresses), "sha256": digest.hexdigest()}
write`);
  const reply = write(() => createReplyBuffer(pyodide.FS)).toJs({ dict_converter: Object.fromEntries });
  // Every chunk sat above 2 GiB, so the test reads the addresses the bug was about.
  assert.ok(reply.chunks > 2);
  assert.equal(reply.low, 0);
  const result = { metadata: reply.buffer, results: [] };
  finishReplyBuffers(result);
  const bytes = new Uint8Array(result.metadata);
  assert.equal(createHash('sha256').update(bytes).digest('hex'), reply.sha256);
  const text = new TextDecoder('utf-8', { fatal: true }).decode(bytes);
  assert.equal(JSON.parse(text).items[0].features[59999].product, 'µ-protein 日本 \u{1F9EC}');
  // The chunk file is gone once its chunk is taken.
  assert.equal(pyodide.FS.analyzePath(reply.buffer.path).exists, false);
} finally {
  release();
}

console.log('reply buffer high-heap tests passed');
