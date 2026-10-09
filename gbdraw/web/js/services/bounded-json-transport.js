// @ts-check
const CLONE_LIMIT = 128 * 1024;
const STRING_UNITS = 128 * 1024;
const BYTE_CHUNK = 256 * 1024;
const PART_OVERHEAD = 64;
const cloneUpperBound = (value, remaining) => {
  if (value instanceof Uint8Array) return remaining + 1;
  if (value === null) return 4;
  if (typeof value === 'string') return 2 + value.length * 6;
  if (typeof value !== 'object') return 25;
  let bytes = 2;
  for (const key of Object.keys(value)) {
    bytes += Array.isArray(value) ? 1 : key.length * 6 + 4;
    if (bytes > remaining) return bytes;
    bytes += cloneUpperBound(value[key], remaining - bytes);
    if (bytes > remaining) return bytes;
  }
  return bytes;
};

// Workers may relinquish an exclusively owned parsed result after each ACK.
// Borrowed values remain unchanged; the receiver always gets the complete graph.
// Consecutive parts share one message up to CLONE_LIMIT (a full byte or string chunk travels
// alone), so messages and ACKs follow the payload size, not its number of keys.
export const sendBoundedJson = async (value, sendPart, { consume = false } = {}) => {
  /** @type {any[]} */
  let parts = [];
  /** @type {Transferable[]} */
  let transfers = [];
  /** @type {(() => void)[]} */
  const releases = [];
  let bytes = 0;
  const flush = async () => {
    if (parts.length) {
      const message = parts.length === 1 ? parts[0] : { kind: 'parts', parts };
      const buffers = transfers;
      parts = [];
      transfers = [];
      bytes = 0;
      await sendPart(message, buffers);
    }
    // Every queued release belongs to a part acknowledged by now.
    releases.splice(0).forEach((release) => release());
  };
  const emit = async (part, size, buffers = []) => {
    if (parts.length && bytes + size > CLONE_LIMIT) await flush();
    parts.push(part);
    transfers.push(...buffers);
    bytes += size;
  };
  const walk = async (current, path) => {
    const pathBytes = PART_OVERHEAD + cloneUpperBound(path, CLONE_LIMIT);
    if (current instanceof Uint8Array) {
      await emit({ kind: 'bytes-start', path, length: current.byteLength }, pathBytes);
      for (let offset = 0; offset < current.byteLength; offset += BYTE_CHUNK) {
        const chunk = current.slice(offset, offset + BYTE_CHUNK);
        await emit({ kind: 'bytes-chunk', offset, bytes: chunk }, PART_OVERHEAD + chunk.byteLength, [chunk.buffer]);
      }
      await emit({ kind: 'bytes-end' }, PART_OVERHEAD);
      return;
    }
    const size = cloneUpperBound(current, CLONE_LIMIT);
    if (size <= CLONE_LIMIT) {
      await emit({ kind: 'value', path, value: current, whole: path.length === 0 }, pathBytes + size);
    } else if (typeof current === 'string') {
      await emit({ kind: 'string-start', path }, pathBytes);
      for (let start = 0; start < current.length; start += STRING_UNITS) {
        const units = new Uint16Array(Math.min(STRING_UNITS, current.length - start));
        for (let index = 0; index < units.length; index++) units[index] = current.charCodeAt(start + index);
        await emit({ kind: 'string-chunk', units }, PART_OVERHEAD + units.byteLength, [units.buffer]);
      }
      await emit({ kind: 'string-end' }, PART_OVERHEAD);
    } else if (Array.isArray(current)) {
      await emit({ kind: 'value', path, value: [] }, pathBytes);
      for (let index = 0; index < current.length;) {
        const batch = [];
        let batchBytes = 2;
        while (index + batch.length < current.length && batch.length < 64) {
          const item = current[index + batch.length];
          const itemBytes = cloneUpperBound(item, CLONE_LIMIT - batchBytes);
          if (batchBytes + itemBytes > CLONE_LIMIT) break;
          batch.push(item);
          batchBytes += itemBytes + 1;
        }
        if (batch.length) {
          const start = index;
          await emit({ kind: 'batch', path, index, value: batch }, pathBytes + batchBytes);
          if (consume) releases.push(() => current.fill(null, start, start + batch.length));
          index += batch.length;
        } else {
          const at = index;
          await walk(current[at], [...path, at]);
          if (consume) releases.push(() => { current[at] = null; });
          index++;
        }
      }
    } else {
      await emit({ kind: 'value', path, value: {} }, pathBytes);
      for (const key of Object.keys(current)) {
        await walk(current[key], [...path, key]);
        if (consume) releases.push(() => { delete current[key]; });
      }
    }
  };
  await walk(value, []);
  await flush();
};

export const createBoundedJsonReceiver = () => {
  let candidate;
  /** @type {{ path: any, chunks: string[] } | null} */
  let pendingString = null;
  /** @type {{ bytes: Uint8Array, offset: number } | null} */
  let pendingBytes = null;
  const resolvePath = (path) => {
    if (!Array.isArray(path)) throw new Error('Invalid bounded JSON transport path.');
    let value = candidate;
    for (const key of path) {
      if (!value || !Object.hasOwn(value, key)) throw new Error('Invalid bounded JSON transport parent.');
      value = value[key];
    }
    return value;
  };
  const assignValue = (path, value) => {
    if (!Array.isArray(path)) throw new Error('Invalid bounded JSON transport path.');
    if (path.length === 0) { candidate = value; return; }
    const parent = resolvePath(path.slice(0, -1));
    // Preserve own unsafe keys for the admission validator without invoking prototype setters.
    Object.defineProperty(parent, path.at(-1), { value, enumerable: true, writable: true, configurable: true });
  };
  const receiveOne = (part) => {
    if (part.kind === 'value') assignValue(part.path, part.value);
    else if (part.kind === 'bytes-start') {
      if (pendingBytes || !Number.isSafeInteger(part.length) || part.length < 0) throw new Error('Invalid bounded bytes.');
      const bytes = new Uint8Array(part.length);
      pendingBytes = { bytes, offset: 0 };
      assignValue(part.path, bytes);
    } else if (part.kind === 'bytes-chunk') {
      if (!pendingBytes || !(part.bytes instanceof Uint8Array)
        || part.bytes.byteLength > BYTE_CHUNK || part.offset !== pendingBytes.offset
        || part.offset + part.bytes.byteLength > pendingBytes.bytes.byteLength) {
        throw new Error('Invalid bounded bytes chunk.');
      }
      pendingBytes.bytes.set(part.bytes, part.offset);
      pendingBytes.offset += part.bytes.byteLength;
    } else if (part.kind === 'bytes-end') {
      if (!pendingBytes || pendingBytes.offset !== pendingBytes.bytes.byteLength) throw new Error('Incomplete bounded bytes.');
      pendingBytes = null;
    }
    else if (part.kind === 'batch') {
      const target = resolvePath(part.path);
      if (!Array.isArray(target) || part.index !== target.length || !Array.isArray(part.value)) {
        throw new Error('Invalid bounded JSON transport batch.');
      }
      target.push(...part.value);
    } else if (part.kind === 'string-start') {
      if (pendingString) throw new Error('Invalid bounded JSON transport string.');
      pendingString = { path: part.path, chunks: [] };
    } else if (part.kind === 'string-chunk') {
      if (!pendingString || !(part.units instanceof Uint16Array)) throw new Error('Invalid bounded JSON transport string.');
      const pieces = [];
      for (let index = 0; index < part.units.length; index += 8192) {
        pieces.push(String.fromCharCode(...part.units.subarray(index, index + 8192)));
      }
      pendingString.chunks.push(pieces.join(''));
    } else if (part.kind === 'string-end') {
      if (!pendingString) throw new Error('Invalid bounded JSON transport string.');
      assignValue(pendingString.path, pendingString.chunks.join(''));
      pendingString = null;
    } else throw new Error('Invalid bounded JSON transport part.');
  };
  // One message carries one part or a `parts` list of them, never a nested list.
  const receivePart = (part) => {
    if (part.kind === 'parts' && Array.isArray(part.parts)) part.parts.forEach(receiveOne);
    else receiveOne(part);
  };
  return { receivePart, getValue: () => {
    if (candidate === undefined || pendingString || pendingBytes) throw new Error('Incomplete bounded JSON transport.');
    return candidate;
  } };
};
