const GZIP_MAGIC = Object.freeze([0x1f, 0x8b]);
const MAX_SESSION_FILE_BYTES = 200 * 1024 * 1024;
const MAX_EXPANDED_SESSION_BYTES = 512 * 1024 * 1024;

const downloadBlob = (blob, filename) => {
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement('a');
  anchor.href = url;
  anchor.download = filename;
  anchor.click();
  URL.revokeObjectURL(url);
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

export const downloadCompressedSession = async (data, filename) => {
  if (typeof CompressionStream !== 'function') {
    throw new Error('This browser does not support gzip session export.');
  }
  const json = JSON.stringify(data);
  const source = new Blob([json], { type: 'application/json' });
  const compressedStream = source.stream().pipeThrough(new CompressionStream('gzip'));
  const compressed = await new Response(compressedStream).blob();
  downloadBlob(compressed.slice(0, compressed.size, 'application/gzip'), filename);
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
