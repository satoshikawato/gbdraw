const GZIP_MAGIC = Object.freeze([0x1f, 0x8b]);
const MAX_SESSION_FILE_BYTES = 200 * 1024 * 1024;
const MAX_EXPANDED_SESSION_BYTES = 512 * 1024 * 1024;
const JSON_CHUNK_TARGET_BYTES = 256 * 1024;
const JSON_CHUNKS_PER_TASK_YIELD = 8;
export const SESSION_DOWNLOAD_CONFIRM_THRESHOLD_BYTES = 50 * 1024 * 1024;

export const confirmLargeSessionBlob = (
  blob,
  confirmFn = globalThis.confirm
) => {
  if (!(blob instanceof Blob)) {
    throw new TypeError('Session download confirmation requires a Blob.');
  }
  if (blob.size <= SESSION_DOWNLOAD_CONFIRM_THRESHOLD_BYTES) return true;
  if (typeof confirmFn !== 'function') return true;
  return Boolean(confirmFn(
    `Compressed session size is ${(blob.size / (1024 * 1024)).toFixed(1)} MB. Continue?`
  ));
};

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
        throw new Error('Expanded session file is too large.');
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

const isGzipFile = async (file) => {
  if (file.size < GZIP_MAGIC.length) return false;
  const header = new Uint8Array(await file.slice(0, GZIP_MAGIC.length).arrayBuffer());
  return GZIP_MAGIC.every((byte, index) => header[index] === byte);
};

const isPlainObject = (value) => {
  if (!value || typeof value !== 'object') return false;
  const prototype = Object.getPrototypeOf(value);
  return prototype === Object.prototype || prototype === null;
};

const unsupportedJsonValue = (value) => {
  const type = typeof value;
  return type === 'undefined' || type === 'function' || type === 'symbol';
};

const boundedJsonUpperBytes = (value, ancestors, container, limit) => {
  if (value === null) return limit >= 4 ? 4 : -1;
  const type = typeof value;
  if (type === 'string') {
    const upperBytes = 2 + (value.length * 6);
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
    || (!toJsonDescriptor && typeof value.toJSON === 'function')
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

function* jsonFragments(value, ancestors = new Set(), container = 'root') {
  if (value === null) {
    yield 'null';
    return;
  }

  const type = typeof value;
  if (type === 'string') {
    yield* quotedJsonFragments(value);
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
  if (typeof value?.toJSON === 'function') {
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

const jsonByteStream = (data) => {
  const fragments = jsonFragments(data);
  const encoder = new TextEncoder();
  let fragment = '';
  let fragmentOffset = 0;
  let complete = false;
  let chunksSinceTaskYield = 0;

  return new ReadableStream({
    async pull(controller) {
      if (chunksSinceTaskYield >= JSON_CHUNKS_PER_TASK_YIELD) {
        chunksSinceTaskYield = 0;
        await new Promise((resolve) => setTimeout(resolve, 0));
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
      if (byteLength > 0) {
        controller.enqueue(chunk.subarray(0, byteLength));
        chunksSinceTaskYield += 1;
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

export const compressSessionData = async (data) => {
  if (typeof CompressionStream !== 'function') {
    throw new Error('This browser does not support gzip session export.');
  }
  const compressedStream = jsonByteStream(data)
    .pipeThrough(new CompressionStream('gzip'));
  return new Response(compressedStream, {
    headers: { 'Content-Type': 'application/gzip' }
  }).blob();
};

export const readSessionText = async (file) => {
  if (file.size > MAX_SESSION_FILE_BYTES) {
    throw new Error('Session file is too large.');
  }
  if (!(await isGzipFile(file))) {
    return readUtf8Stream(file.stream(), MAX_SESSION_FILE_BYTES);
  }
  if (typeof DecompressionStream !== 'function') {
    throw new Error('This browser does not support gzip session import.');
  }
  const expanded = file.stream().pipeThrough(new DecompressionStream('gzip'));
  return readUtf8Stream(expanded, MAX_EXPANDED_SESSION_BYTES);
};
