// Disposable S01 measurement only. The application never imports this Worker.
import { readSessionText } from '/gbdraw/web/js/services/session-file.js';

const LIMIT = 128 * 1024;
const STRING_BYTES = 256 * 1024;
let acknowledge;
let codec;
// Observe the unchanged codec's stream spans; do not reproduce its parser/limits.
const tap = (stream, kind) => {
  const metric = {kind, start: performance.now(), bytes: 0, chunks: 0};
  codec.streams.push(metric);
  return stream.pipeThrough(new TransformStream({
    transform(chunk, controller) {metric.bytes += chunk.byteLength; metric.chunks++; controller.enqueue(chunk);},
    flush() {metric.wallMs = performance.now() - metric.start;}
  }));
};
const blobStream = Blob.prototype.stream;
Blob.prototype.stream = function(...args) {return tap(blobStream.apply(this, args), 'file-read');};
const NativeDecompress = self.DecompressionStream;
self.DecompressionStream = class extends NativeDecompress {
  constructor(...args) {super(...args); return {writable: this.writable, readable: tap(this.readable, 'decompress')};}
};
const decode = TextDecoder.prototype.decode;
TextDecoder.prototype.decode = function(...args) {
  const start = performance.now(); const value = decode.apply(this, args);
  codec.decodeCpuMs += performance.now() - start;
  codec.decodedBytes += args[0]?.byteLength || 0;
  return value;
};
self.onmessage = async ({ data }) => {
  if (data.ack) { acknowledge?.(); return; }
  const { id, file, mode } = data;
  const stages = {};
  codec = {streams: [], decodeCpuMs: 0, decodedBytes: 0};
  let messages = 0;
  let transferredBytes = 0;
  let clonedUpperBytes = mode === 'whole' ? null : 0;
  let postMessageCpuMs = 0;
  let maximumMessageUpperBytes = mode === 'whole' ? null : 0;
  const send = (message, transfers = []) => new Promise(resolve => {
    acknowledge = resolve;
    messages++;
    const start = performance.now();
    self.postMessage({ id, ...message }, transfers);
    postMessageCpuMs += performance.now() - start;
  });
  // A bounded estimate stops before visiting an unbounded subtree. It is not a codec.
  const upper = (value, budget) => {
    if (value === null) return 4;
    if (typeof value === 'string') return 2 + value.length * 6;
    if (typeof value !== 'object') return 25;
    let bytes = 2;
    for (const key of Object.keys(value)) {
      bytes += Array.isArray(value) ? 1 : key.length * 6 + 4;
      if (bytes > budget) return bytes;
      bytes += upper(value[key], budget - bytes);
      if (bytes > budget) return bytes;
    }
    return bytes;
  };
  const emit = async (path, value) => {
    const size = upper(value, LIMIT);
    if (size <= LIMIT) {
      clonedUpperBytes += size;
      maximumMessageUpperBytes = Math.max(maximumMessageUpperBytes, size);
      await send({ kind: 'value', path, value });
    } else if (typeof value === 'string') {
      await send({ kind: 'string-start', path });
      // Preserve every JSON string code unit, including escaped lone surrogates.
      for (let start = 0; start < value.length; start += STRING_BYTES / 2) {
        const end = Math.min(value.length, start + STRING_BYTES / 2);
        const units = new Uint16Array(end - start);
        for (let index = start; index < end; index++) units[index - start] = value.charCodeAt(index);
        transferredBytes += units.byteLength;
        maximumMessageUpperBytes = Math.max(maximumMessageUpperBytes, units.byteLength);
        await send({ kind: 'string-chunk', path, units }, [units.buffer]);
        if (units.byteLength !== 0) throw new Error('Transfer buffer did not detach');
      }
      await send({ kind: 'string-end', path });
    } else if (Array.isArray(value)) {
      await send({ kind: 'value', path, value: [] });
      for (let index = 0; index < value.length;) {
        const batch = [];
        let batchBytes = 2;
        while (index + batch.length < value.length && batch.length < 64) {
          const entry = value[index + batch.length];
          const bytes = upper(entry, LIMIT - batchBytes);
          if (batchBytes + bytes > LIMIT) break;
          batch.push(entry);
          batchBytes += bytes + 1;
        }
        if (batch.length) {
          clonedUpperBytes += batchBytes;
          maximumMessageUpperBytes = Math.max(maximumMessageUpperBytes, batchBytes);
          await send({ kind: 'batch', path, index, value: batch });
          index += batch.length;
        } else {
          await emit([...path, index], value[index]);
          index++;
        }
      }
    } else {
      await send({ kind: 'value', path, value: {} });
      for (const [key, entry] of Object.entries(value)) await emit([...path, key], entry);
    }
  };
  try {
    const readStart = performance.now();
    const text = await readSessionText(file);
    stages.readDecodeDecompressMs = performance.now() - readStart;
    const parseStart = performance.now();
    const document = JSON.parse(text);
    stages.parseMs = performance.now() - parseStart;
    self.postMessage({ id, kind: 'parsed', stages, expandedCharacters: text.length });
    // Let the independent CDP sampler observe the parsed graph before transfer.
    await new Promise(resolve => setTimeout(resolve, 120));
    const transferStart = performance.now();
    if (mode === 'whole') {
      const postStart = performance.now();
      await send({ kind: 'value', path: [], value: document });
      stages.wholeReplyRoundTripMs = performance.now() - postStart;
    } else if (mode === 'bounded') {
      await emit([], document);
    } else {
      throw new Error('Unknown disposable transport');
    }
    stages.transferAndAssemblyMs = performance.now() - transferStart;
    self.postMessage({ id, kind: 'done', stages, codec, messages, transferredBytes,
                       clonedUpperBytes, maximumMessageUpperBytes, postMessageCpuMs });
  } catch (error) {
    self.postMessage({ id, kind: 'error', stage: 'probe', message: error.message });
  }
};
