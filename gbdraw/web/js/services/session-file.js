// @ts-check
import { diagnosticError } from '../utils/error-normalization.js';

const GZIP_MAGIC = Object.freeze([0x1f, 0x8b]);
const MAX_SESSION_FILE_BYTES = 200 * 1024 * 1024;
// The most JSON a Session may expand to: Load stops reading past it, and Save
// refuses to write more (OV-2240), so a saved file always loads again.
export const MAX_EXPANDED_SESSION_BYTES = 512 * 1024 * 1024;
const LARGEST_INPUT_COUNT = 3;
const LARGEST_INPUT_NAME_CHARS = 80;
const JSON_CHUNK_TARGET_BYTES = 256 * 1024;
const JSON_TASK_BUDGET_MS = 16;
export const SESSION_DOWNLOAD_CONFIRM_THRESHOLD_BYTES = 50 * 1024 * 1024;

/**
 * A compressed Session above the threshold asks before its download;
 * the caller asks with the app's choice dialog.
 * @param {Blob} blob The compressed Session.
 * @returns {number | null} Its size in MB when the download asks first, else null.
 */
export const largeSessionDownloadMegabytes = (blob) => {
  if (!(blob instanceof Blob)) {
    throw new TypeError('Session download confirmation requires a Blob.');
  }
  return blob.size > SESSION_DOWNLOAD_CONFIRM_THRESHOLD_BYTES ? blob.size / (1024 * 1024) : null;
};

/**
 * @param {ReadableStream<Uint8Array>} stream
 * @param {number} maxBytes
 * @returns {Promise<string>}
 */
const readUtf8Stream = async (stream, maxBytes) => {
  const reader = stream.getReader();
  const decoder = new TextDecoder('utf-8', { fatal: true });
  let byteCount = 0;
  let text = '';

  try {
    while (true) {
      const { value, done } = await reader.read();
      if (done) break;
      byteCount += value.byteLength;
      if (byteCount > maxBytes) {
        throw diagnosticError('SESSION_SIZE_LIMIT', {}, { stage: 'read' });
      }
      text += decoder.decode(value, { stream: true });
    }
    return text + decoder.decode();
  } catch (error) {
    await reader.cancel(error).catch(() => {});
    throw error;
  } finally {
    reader.releaseLock();
  }
};

/** @param {Blob} file */
const isGzipFile = async (file) => {
  if (file.size < GZIP_MAGIC.length) return false;
  const header = new Uint8Array(await file.slice(0, GZIP_MAGIC.length).arrayBuffer());
  return GZIP_MAGIC.every((byte, index) => header[index] === byte);
};

/**
 * @param {unknown} value
 * @returns {value is Record<string, unknown>}
 */
const isPlainObject = (value) => {
  if (!value || typeof value !== 'object') return false;
  const prototype = Object.getPrototypeOf(value);
  return prototype === Object.prototype || prototype === null;
};

/** @param {unknown} value */
const unsupportedJsonValue = (value) => {
  const type = typeof value;
  return type === 'undefined' || type === 'function' || type === 'symbol';
};

/** @typedef {'root' | 'array' | 'object'} JsonContainer Where a value sits in the document being sized or written. */

/**
 * An upper bound of the JSON bytes `value` needs, or -1 when it exceeds `limit`
 * or JSON.stringify would not treat it as plain data.
 * @param {unknown} value
 * @param {Set<unknown>} ancestors
 * @param {JsonContainer} container
 * @param {number} limit
 * @returns {number}
 */
const boundedJsonUpperBytes = (value, ancestors, container, limit) => {
  if (value === null) return limit >= 4 ? 4 : -1;
  const type = typeof value;
  if (type === 'string') {
    const upperBytes = 2 + (/** @type {string} */ (value).length * 6); // `type === 'string'` above
    return upperBytes <= limit ? upperBytes : -1;
  }
  if (type === 'boolean') {
    const upperBytes = value ? 4 : 5;
    return upperBytes <= limit ? upperBytes : -1;
  }
  if (type === 'number') return limit >= 25 ? 25 : -1;
  if (type === 'bigint') return -1;
  if (unsupportedJsonValue(value)) {
    if (container !== 'array') return 0;
    return limit >= 4 ? 4 : -1;
  }
  if (!Array.isArray(value) && !isPlainObject(value)) return -1;
  if (ancestors.has(value)) return -1;

  const toJsonDescriptor = Object.getOwnPropertyDescriptor(value, 'toJSON');
  if (
    toJsonDescriptor?.get
    || typeof toJsonDescriptor?.value === 'function'
    || (!toJsonDescriptor && typeof /** @type {{ toJSON?: unknown }} */ (value).toJSON === 'function') // probe of an array or a plain object
  ) return -1;

  ancestors.add(value);
  try {
    let upperBytes = 2;
    if (Array.isArray(value)) {
      for (let index = 0; index < value.length; index += 1) {
        if (index > 0) upperBytes += 1;
        const descriptor = Object.getOwnPropertyDescriptor(value, String(index));
        if (!descriptor && index in value) return -1;
        if (descriptor?.get || descriptor?.set) return -1;
        const itemBytes = boundedJsonUpperBytes(
          descriptor?.value,
          ancestors,
          'array',
          limit - upperBytes
        );
        if (itemBytes < 0) return -1;
        upperBytes += itemBytes;
        if (upperBytes > limit) return -1;
      }
      return upperBytes;
    }

    let emitted = 0;
    for (const key of Object.keys(value)) {
      const descriptor = Object.getOwnPropertyDescriptor(value, key);
      if (!descriptor || descriptor.get || descriptor.set) return -1;
      const entryBytes = boundedJsonUpperBytes(
        descriptor.value,
        ancestors,
        'object',
        limit - upperBytes
      );
      if (entryBytes < 0) return -1;
      if (entryBytes === 0) continue;
      upperBytes += (emitted > 0 ? 1 : 0) + 3 + (key.length * 6) + entryBytes;
      if (upperBytes > limit) return -1;
      emitted += 1;
    }
    return upperBytes;
  } finally {
    ancestors.delete(value);
  }
};

/**
 * @param {unknown} value
 * @param {Set<unknown>} ancestors
 * @param {JsonContainer} container
 * @returns {string | null}
 */
const boundedNativeJson = (value, ancestors, container) => (
  boundedJsonUpperBytes(
    value,
    ancestors,
    container,
    Math.floor(JSON_CHUNK_TARGET_BYTES / 2)
  ) >= 0
    ? JSON.stringify(value)
    : null
);

/**
 * @param {string} value
 * @returns {Generator<string, void, undefined>}
 */
function* quotedJsonFragments(value) {
  yield '"';
  const sliceTarget = Math.max(1, Math.floor(JSON_CHUNK_TARGET_BYTES / 8));
  for (let start = 0; start < value.length;) {
    let end = Math.min(value.length, start + sliceTarget);
    if (
      end < value.length
      && value.charCodeAt(end - 1) >= 0xd800
      && value.charCodeAt(end - 1) <= 0xdbff
      && value.charCodeAt(end) >= 0xdc00
      && value.charCodeAt(end) <= 0xdfff
    ) {
      end += 1;
    }
    const encoded = JSON.stringify(value.slice(start, end));
    yield encoded.slice(1, -1);
    start = end;
  }
  yield '"';
}

/**
 * @param {unknown} value
 * @param {Set<unknown>} [ancestors]
 * @param {JsonContainer} [container]
 * @returns {Generator<string, void, undefined>}
 */
function* jsonFragments(value, ancestors = new Set(), container = 'root') {
  if (value === null) {
    yield 'null';
    return;
  }

  const type = typeof value;
  if (type === 'string') {
    yield* quotedJsonFragments(/** @type {string} */ (value)); // `type === 'string'` above
    return;
  }
  if (type === 'boolean') {
    yield value ? 'true' : 'false';
    return;
  }
  if (type === 'number') {
    yield Number.isFinite(value) ? (Object.is(value, -0) ? '0' : String(value)) : 'null';
    return;
  }
  if (type === 'bigint') {
    throw new TypeError('BigInt values cannot be written to a Session file.');
  }
  if (unsupportedJsonValue(value)) {
    if (container === 'array') yield 'null';
    else if (container === 'root') {
      throw new TypeError('Session data must be a JSON-compatible value.');
    }
    return;
  }
  if (typeof /** @type {{ toJSON?: unknown } | undefined} */ (value)?.toJSON === 'function') {
    throw new TypeError('Session data contains an unsupported object value.');
  }
  if (!Array.isArray(value) && !isPlainObject(value)) {
    throw new TypeError('Session data contains an unsupported object value.');
  }
  if (ancestors.has(value)) {
    throw new TypeError('Session data contains a cyclic value.');
  }

  const nativeJson = boundedNativeJson(value, ancestors, container);
  if (nativeJson !== null) {
    yield nativeJson;
    return;
  }

  ancestors.add(value);
  try {
    if (Array.isArray(value)) {
      yield '[';
      for (let index = 0; index < value.length; index += 1) {
        if (index > 0) yield ',';
        yield* jsonFragments(value[index], ancestors, 'array');
      }
      yield ']';
      return;
    }

    yield '{';
    let emitted = 0;
    for (const key of Object.keys(value)) {
      const entry = value[key];
      if (unsupportedJsonValue(entry)) continue;
      if (emitted > 0) yield ',';
      yield* quotedJsonFragments(key);
      yield ':';
      yield* jsonFragments(entry, ancestors, 'object');
      emitted += 1;
    }
    yield '}';
  } finally {
    ancestors.delete(value);
  }
}

/**
 * The file names of the Session's largest inputs, largest first.
 * @param {unknown} resources The Session's `resources` table.
 * @returns {string}
 */
const largestInputNames = (resources) => {
  if (!isPlainObject(resources)) return '';
  /** @type {{ name: string, size: number }[]} */
  const inputs = [];
  for (const entry of Object.values(resources)) {
    if (!isPlainObject(entry) || typeof entry.name !== 'string' || !entry.name) continue;
    inputs.push({
      name: entry.name.slice(0, LARGEST_INPUT_NAME_CHARS),
      size: Number.isSafeInteger(entry.size) ? /** @type {number} */ (entry.size) : 0
    });
  }
  return inputs.sort((a, b) => b.size - a.size).slice(0, LARGEST_INPUT_COUNT)
    .map((input) => input.name).join(', ');
};

/**
 * @param {Record<string, unknown>} data The Session document.
 * @param {number} sessionBytes Its expanded JSON size.
 * @param {number} sessionByteLimit The cap it exceeds.
 */
const sessionSaveSizeLimit = (data, sessionBytes, sessionByteLimit) => diagnosticError(
  'SESSION_SIZE_LIMIT',
  { sessionBytes, sessionByteLimit, largestInputs: largestInputNames(data.resources) },
  { stage: 'write', operation: 'session-save' }
);

/**
 * The Session's JSON as bytes. Past `maxBytes` it writes nothing more, counts
 * the rest, and fails with the whole size.
 * @param {Record<string, unknown>} data
 * @param {number} maxBytes
 * @returns {ReadableStream<Uint8Array<ArrayBuffer>>}
 */
const jsonByteStream = (data, maxBytes) => {
  const fragments = jsonFragments(data);
  const encoder = new TextEncoder();
  let fragment = '';
  let fragmentOffset = 0;
  let complete = false;
  let writtenBytes = 0;
  let taskDeadline = performance.now() + JSON_TASK_BUDGET_MS;

  return new ReadableStream({
    async pull(controller) {
      if (performance.now() >= taskDeadline) {
        await new Promise((resolve) => setTimeout(resolve, 0));
        taskDeadline = performance.now() + JSON_TASK_BUDGET_MS;
      }
      const chunk = new Uint8Array(JSON_CHUNK_TARGET_BYTES);
      let byteLength = 0;
      while (byteLength < JSON_CHUNK_TARGET_BYTES) {
        if (fragmentOffset >= fragment.length) {
          const next = fragments.next();
          if (next.done) {
            complete = true;
            break;
          }
          fragment = next.value;
          fragmentOffset = 0;
        }
        const encoded = encoder.encodeInto(
          fragment.slice(fragmentOffset),
          chunk.subarray(byteLength)
        );
        if (encoded.read === 0 && encoded.written === 0) {
          break;
        }
        fragmentOffset += encoded.read;
        byteLength += encoded.written;
      }
      writtenBytes += byteLength;
      if (writtenBytes > maxBytes) {
        let totalBytes = writtenBytes + encoder.encode(fragment.slice(fragmentOffset)).byteLength;
        for (let next = fragments.next(); !next.done; next = fragments.next()) {
          totalBytes += encoder.encode(next.value).byteLength;
          if (performance.now() >= taskDeadline) {
            await new Promise((resolve) => setTimeout(resolve, 0));
            taskDeadline = performance.now() + JSON_TASK_BUDGET_MS;
          }
        }
        complete = true;
        controller.error(sessionSaveSizeLimit(data, totalBytes, maxBytes));
        return;
      }
      if (byteLength > 0) {
        controller.enqueue(chunk.subarray(0, byteLength));
      }
      if (complete) controller.close();
    },
    cancel() {
      fragments.return?.();
      fragment = '';
      fragmentOffset = 0;
      complete = true;
    }
  });
};

/**
 * @param {Record<string, any>} data The Session document.
 * @param {{ maxExpandedBytes?: number }} [options] `maxExpandedBytes` stands in for the cap in tests.
 * @returns {Promise<Blob>} The gzip Session; it rejects with `SESSION_SIZE_LIMIT` (stage `write`)
 *   when the JSON is over the cap.
 */
export const compressSessionData = async (data, { maxExpandedBytes = MAX_EXPANDED_SESSION_BYTES } = {}) => {
  if (typeof CompressionStream !== 'function') {
    throw diagnosticError('SESSION_BROWSER_UNSUPPORTED', {}, { stage: 'transport' });
  }
  const compressedStream = jsonByteStream(data, maxExpandedBytes)
    .pipeThrough(new CompressionStream('gzip'));
  return new Response(compressedStream, {
    headers: { 'Content-Type': 'application/gzip' }
  }).blob();
};

/**
 * @param {Blob} file A plain or gzip-compressed Session file.
 * @returns {Promise<string>} The Session JSON text; the caller parses and validates it.
 */
export const readSessionText = async (file) => {
  if (file.size > MAX_SESSION_FILE_BYTES) {
    throw diagnosticError('SESSION_SIZE_LIMIT', {}, { stage: 'read' });
  }
  if (!(await isGzipFile(file))) {
    return readUtf8Stream(file.stream(), MAX_SESSION_FILE_BYTES);
  }
  if (typeof DecompressionStream !== 'function') {
    throw diagnosticError('SESSION_BROWSER_UNSUPPORTED', {}, { stage: 'read' });
  }
  const expanded = file.stream().pipeThrough(new DecompressionStream('gzip'));
  return readUtf8Stream(expanded, MAX_EXPANDED_SESSION_BYTES);
};
