// @ts-check
const CLONE_LIMIT = 128 * 1024;
const STRING_UNITS = 128 * 1024;
const BYTE_CHUNK = 256 * 1024;
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
export const sendBoundedJson = async (value, sendPart, path = [], { consume = false } = {}) => {
  if (value instanceof Uint8Array) {
    await sendPart({ kind: 'bytes-start', path, length: value.byteLength });
    for (let offset = 0; offset < value.byteLength; offset += BYTE_CHUNK) {
      const bytes = value.slice(offset, offset + BYTE_CHUNK);
      await sendPart({ kind: 'bytes-chunk', offset, bytes }, [bytes.buffer]);
    }
    await sendPart({ kind: 'bytes-end' });
  } else if (cloneUpperBound(value, CLONE_LIMIT) <= CLONE_LIMIT) {
    await sendPart({ kind: 'value', path, value, whole: path.length === 0 });
  } else if (typeof value === 'string') {
    await sendPart({ kind: 'string-start', path });
    for (let start = 0; start < value.length; start += STRING_UNITS) {
      const units = new Uint16Array(Math.min(STRING_UNITS, value.length - start));
      for (let index = 0; index < units.length; index++) units[index] = value.charCodeAt(start + index);
      await sendPart({ kind: 'string-chunk', units }, [units.buffer]);
    }
    await sendPart({ kind: 'string-end' });
  } else if (Array.isArray(value)) {
    await sendPart({ kind: 'value', path, value: [] });
    for (let index = 0; index < value.length;) {
      const batch = [];
      let bytes = 2;
      while (index + batch.length < value.length && batch.length < 64) {
        const item = value[index + batch.length];
        const size = cloneUpperBound(item, CLONE_LIMIT - bytes);
        if (bytes + size > CLONE_LIMIT) break;
        batch.push(item);
        bytes += size + 1;
      }
      if (batch.length) {
        await sendPart({ kind: 'batch', path, index, value: batch });
        if (consume) value.fill(null, index, index + batch.length);
        index += batch.length;
      } else {
        await sendBoundedJson(value[index], sendPart, [...path, index], { consume });
        if (consume) value[index] = null;
        index++;
      }
    }
  } else {
    await sendPart({ kind: 'value', path, value: {} });
    for (const key of Object.keys(value)) {
      await sendBoundedJson(value[key], sendPart, [...path, key], { consume });
      if (consume) delete value[key];
    }
  }
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
  const receivePart = (part) => {
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
  return { receivePart, getValue: () => {
    if (candidate === undefined || pendingString || pendingBytes) throw new Error('Incomplete bounded JSON transport.');
    return candidate;
  } };
};
