// @ts-check
const CLONE_LIMIT = 128 * 1024;
const STRING_UNITS = 128 * 1024;
const BYTE_CHUNK = 256 * 1024;
const PART_OVERHEAD = 64;
const UTF8 = new TextEncoder();

/** @typedef {null | boolean | number | string | Uint8Array | BoundedJsonValue[] | { [key: string]: BoundedJsonValue }} BoundedJsonValue */
/** @typedef {(string | number)[]} BoundedJsonPath */
/** @typedef {{ [key: string]: BoundedJsonValue }} BoundedJsonRecord An object, or an array read through its own keys. */
/**
 * One transport message. The receiver validates every field, so each member lists what the sender writes.
 * @typedef {{ kind: 'value', path: BoundedJsonPath, value: BoundedJsonValue, whole?: boolean }
 *   | { kind: 'bytes-start', path: BoundedJsonPath, length: number }
 *   | { kind: 'bytes-chunk', offset: number, bytes: Uint8Array }
 *   | { kind: 'bytes-end' }
 *   | { kind: 'batch', path: BoundedJsonPath, index: number, value: BoundedJsonValue[] }
 *   | { kind: 'string-start', path: BoundedJsonPath, encoding: BoundedStringEncoding }
 *   | { kind: 'string-chunk', units: Uint16Array }
 *   | { kind: 'string-bytes', bytes: Uint8Array }
 *   | { kind: 'string-end' }} BoundedJsonSinglePart
 */
/** @typedef {BoundedJsonSinglePart | { kind: 'parts', parts: BoundedJsonSinglePart[] }} BoundedJsonPart */
/**
 * How a long string travels: a well-formed string as UTF-8 `string-bytes` chunks, a string with a lone surrogate
 * (which UTF-8 cannot carry) as UTF-16 `string-chunk` code units.
 * @typedef {'utf-8' | 'utf-16'} BoundedStringEncoding
 */

/**
 * @param {BoundedJsonValue} value
 * @param {number} remaining
 * @returns {number}
 */
const cloneUpperBound = (value, remaining) => {
  if (value instanceof Uint8Array) return remaining + 1;
  if (value === null) return 4;
  if (typeof value === 'string') return 2 + value.length * 6;
  if (typeof value !== 'object') return 25;
  let bytes = 2;
  for (const key of Object.keys(value)) {
    bytes += Array.isArray(value) ? 1 : key.length * 6 + 4;
    if (bytes > remaining) return bytes;
    // An array is walked through its own keys, so it is indexed like a record here.
    bytes += cloneUpperBound(/** @type {BoundedJsonRecord} */ (value)[key], remaining - bytes);
    if (bytes > remaining) return bytes;
  }
  return bytes;
};

// Workers may relinquish an exclusively owned parsed result after each ACK.
// Borrowed values remain unchanged; the receiver always gets the complete graph.
// Consecutive parts share one message up to CLONE_LIMIT (a full byte or string chunk travels
// alone), so messages and ACKs follow the payload size, not its number of keys.
/**
 * @param {BoundedJsonValue} value
 * @param {(message: BoundedJsonPart, buffers: Transferable[]) => void | Promise<void>} sendPart
 * @param {{ consume?: boolean }} [options]
 * @returns {Promise<void>}
 */
export const sendBoundedJson = async (value, sendPart, { consume = false } = {}) => {
  /** @type {BoundedJsonSinglePart[]} */
  let parts = [];
  /** @type {Transferable[]} */
  let transfers = [];
  /** @type {(() => void)[]} */
  const releases = [];
  let bytes = 0;
  const flush = async () => {
    if (parts.length) {
      /** @type {BoundedJsonPart} */
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
  /**
   * @param {BoundedJsonSinglePart} part
   * @param {number} size
   * @param {Transferable[]} [buffers]
   */
  const emit = async (part, size, buffers = []) => {
    if (parts.length && bytes + size > CLONE_LIMIT) await flush();
    parts.push(part);
    transfers.push(...buffers);
    bytes += size;
  };
  /**
   * @param {BoundedJsonValue} current
   * @param {BoundedJsonPath} path
   * @returns {Promise<void>}
   */
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
    } else if (typeof current === 'string' && typeof current.isWellFormed === 'function' && current.isWellFormed()) {
      await emit({ kind: 'string-start', path, encoding: 'utf-8' }, pathBytes);
      for (let start = 0; start < current.length;) {
        // A code unit takes at most 3 bytes, so the last chunk holds the rest. `encodeInto` writes whole characters
        // only; a source slice that ends inside a surrogate pair cannot reach its end, since every unit before it
        // takes at least one of the chunk's bytes.
        const bytes = new Uint8Array(Math.min(BYTE_CHUNK, (current.length - start) * 3));
        const { read, written } = UTF8.encodeInto(current.slice(start, start + bytes.length), bytes);
        start += read;
        await emit({ kind: 'string-bytes', bytes: bytes.subarray(0, written) }, PART_OVERHEAD + written, [bytes.buffer]);
      }
      await emit({ kind: 'string-end' }, PART_OVERHEAD);
    } else if (typeof current === 'string') {
      await emit({ kind: 'string-start', path, encoding: 'utf-16' }, pathBytes);
      for (let start = 0; start < current.length; start += STRING_UNITS) {
        const units = new Uint16Array(Math.min(STRING_UNITS, current.length - start));
        for (let index = 0; index < units.length; index++) units[index] = current.charCodeAt(start + index);
        await emit({ kind: 'string-chunk', units }, PART_OVERHEAD + units.byteLength, [units.buffer]);
      }
      await emit({ kind: 'string-end' }, PART_OVERHEAD);
    } else if (Array.isArray(current)) {
      await emit({ kind: 'value', path, value: [] }, pathBytes);
      for (let index = 0; index < current.length;) {
        /** @type {BoundedJsonValue[]} */
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
      // A scalar fits CLONE_LIMIT and took the first branch, so only an object reaches this one.
      for (const key of Object.keys(/** @type {BoundedJsonRecord} */ (current))) {
        await walk(/** @type {BoundedJsonRecord} */ (current)[key], [...path, key]);
        if (consume) releases.push(() => { delete (/** @type {BoundedJsonRecord} */ (current))[key]; });
      }
    }
  };
  await walk(value, []);
  await flush();
};

/**
 * Decodes one chunk of a UTF-8 string; without `bytes`, ends the string. A chunk may end inside a character.
 * @param {TextDecoder} decoder
 * @param {Uint8Array} [bytes]
 * @returns {string | null} null when the bytes are malformed
 */
const decodeStringBytes = (decoder, bytes) => {
  try {
    return bytes ? decoder.decode(bytes, { stream: true }) : decoder.decode();
  } catch {
    return null;
  }
};

/**
 * @returns {{ receivePart: (part: BoundedJsonPart) => void, getValue: () => BoundedJsonValue }}
 */
export const createBoundedJsonReceiver = () => {
  /** @type {BoundedJsonValue | undefined} */
  let candidate;
  /** @type {{ path: BoundedJsonPath, chunks: string[], decoder: TextDecoder | null } | null} */
  let pendingString = null;
  /** @type {{ bytes: Uint8Array, offset: number } | null} */
  let pendingBytes = null;
  /** @param {BoundedJsonPath} path */
  const resolvePath = (path) => {
    if (!Array.isArray(path)) throw new Error('Invalid bounded JSON transport path.');
    // A path addresses containers. A primitive parent has no own key, or fails when `assignValue` defines on it,
    // and rejects the part, so the walk treats every value as a container.
    let value = /** @type {BoundedJsonRecord | undefined} */ (candidate);
    for (const key of path) {
      if (!value || !Object.hasOwn(value, key)) throw new Error('Invalid bounded JSON transport parent.');
      value = /** @type {BoundedJsonRecord} */ (value[key]);
    }
    return value;
  };
  /**
   * @param {BoundedJsonPath} path
   * @param {BoundedJsonValue} value
   */
  const assignValue = (path, value) => {
    if (!Array.isArray(path)) throw new Error('Invalid bounded JSON transport path.');
    if (path.length === 0) { candidate = value; return; }
    const parent = resolvePath(path.slice(0, -1));
    // Preserve own unsafe keys for the admission validator without invoking prototype setters.
    // `path` is not empty here; `defineProperty` throws for a missing or primitive parent, which rejects the part.
    Object.defineProperty(/** @type {BoundedJsonRecord} */ (parent), /** @type {string | number} */ (path.at(-1)), { value, enumerable: true, writable: true, configurable: true });
  };
  /** @param {BoundedJsonPart} part A nested `parts` message reaches the final `else` and is rejected. */
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
      if (pendingString || (part.encoding !== 'utf-8' && part.encoding !== 'utf-16')) {
        throw new Error('Invalid bounded JSON transport string.');
      }
      // `ignoreBOM` keeps a leading U+FEFF, which is part of the string; `fatal` rejects malformed bytes.
      const decoder = part.encoding === 'utf-8' ? new TextDecoder('utf-8', { fatal: true, ignoreBOM: true }) : null;
      pendingString = { path: part.path, chunks: [], decoder };
    } else if (part.kind === 'string-bytes' || part.kind === 'string-chunk') {
      /** @type {string | null} */
      let piece = null;
      if (part.kind === 'string-bytes') {
        if (pendingString && pendingString.decoder && part.bytes instanceof Uint8Array && part.bytes.byteLength <= BYTE_CHUNK) {
          piece = decodeStringBytes(pendingString.decoder, part.bytes);
        }
      } else if (pendingString && !pendingString.decoder && part.units instanceof Uint16Array) {
        const pieces = [];
        for (let index = 0; index < part.units.length; index += 8192) {
          pieces.push(String.fromCharCode(...part.units.subarray(index, index + 8192)));
        }
        piece = pieces.join('');
      }
      if (!pendingString || piece === null) throw new Error('Invalid bounded JSON transport string.');
      pendingString.chunks.push(piece);
    } else if (part.kind === 'string-end') {
      // The final decode rejects a character left incomplete by the last chunk.
      const tail = pendingString && pendingString.decoder ? decodeStringBytes(pendingString.decoder) : '';
      if (!pendingString || tail === null) throw new Error('Invalid bounded JSON transport string.');
      pendingString.chunks.push(tail);
      assignValue(pendingString.path, pendingString.chunks.join(''));
      pendingString = null;
    } else throw new Error('Invalid bounded JSON transport part.');
  };
  // One message carries one part or a `parts` list of them, never a nested list.
  /** @param {BoundedJsonPart} part */
  const receivePart = (part) => {
    if (part.kind === 'parts' && Array.isArray(part.parts)) part.parts.forEach(receiveOne);
    else receiveOne(part);
  };
  return { receivePart, getValue: () => {
    if (candidate === undefined || pendingString || pendingBytes) throw new Error('Incomplete bounded JSON transport.');
    return candidate;
  } };
};
