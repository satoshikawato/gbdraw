import {
  isSessionResourceFileView,
  readSessionResourceBytes,
  readSessionResourceText,
  sessionResourceSource
} from './session-resource-backing.js';
import {
  asBytes,
  base64ToBytes,
  bytesToBase64,
  bytesToText,
  sha256Hex,
  textToBase64,
  textToBytes
} from './byte-utils.js';

export {
  base64ToBytes,
  bytesToBase64,
  bytesToText,
  sha256Hex,
  textToBase64,
  textToBytes
};

const fileByteCache = new WeakMap();

const asCacheKey = (file) => (
  file && (typeof file === 'object' || typeof file === 'function') ? file : null
);

export const readFileBytes = (file) => {
  try {
    const key = asCacheKey(file);
    const lazyBytes = key && isSessionResourceFileView(file)
      ? readSessionResourceBytes(file)
      : null;
    if (!key || (!lazyBytes && typeof file.arrayBuffer !== 'function')) {
      return Promise.reject(new TypeError(
        'A File-like object with arrayBuffer() is required.'
      ));
    }
    let pending = fileByteCache.get(key);
    if (!pending) {
      pending = lazyBytes || Promise.resolve()
        .then(() => file.arrayBuffer())
        .then(asBytes);
      fileByteCache.set(key, pending);
      if (!lazyBytes) {
        pending.catch(() => {
          if (fileByteCache.get(key) === pending) fileByteCache.delete(key);
        });
      }
    }
    return pending;
  } catch (error) {
    return Promise.reject(error);
  }
};

export const readFileText = async (file) => {
  const lazyText = isSessionResourceFileView(file)
    ? readSessionResourceText(file)
    : null;
  if (lazyText) return lazyText;
  if (typeof file?.arrayBuffer === 'function') {
    return bytesToText(await readFileBytes(file));
  }
  if (typeof file?.text === 'function') return file.text();
  throw new TypeError('A File-like object with arrayBuffer() or text() is required.');
};

export const getSessionResourceSource = (file) => sessionResourceSource(file);

export const cloneFileBytesForTransfer = async (file) => {
  const bytes = await readFileBytes(file);
  return bytes.slice().buffer;
};
