// A Generate reply crosses from Python to the Worker whatever its address in the
// Pyodide heap (OV-2800). Pyodide 0.29.0 reads a Python buffer's address as a
// signed 32-bit number, so a chunk above 2 GiB failed with a RangeError when
// the Worker read it through a buffer view; it now crosses through a file.
// The real vendored Pyodide runs here; the heap is pushed past 2 GiB with blocks
// that are allocated but never touched, so the process stays small.
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

globalThis.self = {};

const root = fileURLToPath(new URL('../../', import.meta.url));
const indexURL = `${root}gbdraw/web/vendor/pyodide/v0.29.0/full/`;
const { loadPyodide } = await import(new URL('../../gbdraw/web/vendor/pyodide/v0.29.0/full/pyodide.mjs', import.meta.url).href);
const { createReplyBuffer, finishReplyBuffers } = await import('../../gbdraw/web/js/workers/diagram-generation-worker.js');

const pyodide = await loadPyodide({ indexURL });
pyodide.FS.writeFile('/home/pyodide/reply_buffers.py', readFileSync(`${root}gbdraw/web_support/reply_buffers.py`));

const HIGH = 2 ** 31;
const highAddress = pyodide.runPython(`
def high_address():
    return id(b"x" * (1 << 20))
high_address`);
// Untouched blocks of about 1 GiB each until a new 1 MiB bytes object lands above 2 GiB.
const pads = [];
while (highAddress() < HIGH) {
  assert.ok(pads.length < 4, 'the Pyodide heap did not grow past 2 GiB');
  const pad = pyodide._module._PyMem_RawMalloc(2 ** 30 - 2 ** 20);
  assert.notEqual(pad, 0);
  pads.push(pad);
}

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
  pads.forEach((pad) => pyodide._module._PyMem_RawFree(pad));
}

console.log('reply buffer high-heap tests passed');
