// @ts-check
import {
  isSessionResourceFileView,
  releaseSessionResourceContent,
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

/** @import { FileLike } from './session-resource-backing.js' */

/**
 * The sha256 hex of a resource's bytes: the Worker's cache key for them, so
 * equal bytes share one token whatever resource IDs name them.
 * @typedef {string & { readonly __brand: 'ResourceContentToken' }} ResourceContentToken
 */

/** @type {WeakMap<object, Promise<Uint8Array>>} */
const fileByteCache = new WeakMap();
// A File's bytes never change, so its token outlives the release of its bytes.
/** @type {WeakMap<object, Promise<ResourceContentToken>>} */
const fileContentTokenCache = new WeakMap();

/**
 * @param {unknown} file
 * @returns {object | null}
 */
const asCacheKey = (file) => (
  file && (typeof file === 'object' || typeof file === 'function') ? file : null
);

/**
 * @param {unknown} file A `File`, a Session resource view, or any value (rejected unless it has `arrayBuffer()`).
 * @returns {Promise<Uint8Array>}
 */
export const readFileBytes = (file) => {
  try {
    const key = asCacheKey(file);
    // A Session resource view caches its bytes in the backing it shares with
    // every other view of the same resource, so a transfer through one view
    // releases them for all (OV-304); the cache below is for ordinary Files.
    // A view has a backing, so its read is never null.
    if (key && isSessionResourceFileView(file)) return /** @type {Promise<Uint8Array>} */ (readSessionResourceBytes(file));
    // `key` is `file` when it is set.
    if (!key || typeof /** @type {{ readonly arrayBuffer?: unknown }} */ (key).arrayBuffer !== 'function') {
      return Promise.reject(new TypeError(
        'A File-like object with arrayBuffer() is required.'
      ));
    }
    let pending = fileByteCache.get(key);
    if (!pending) {
      pending = Promise.resolve()
        // The `typeof` check above holds here; narrowing does not reach the closure.
        .then(() => /** @type {{ arrayBuffer: () => Promise<ArrayBuffer> }} */ (file).arrayBuffer())
        .then(asBytes);
      fileByteCache.set(key, pending);
      pending.catch(() => {
        if (fileByteCache.get(key) === pending) fileByteCache.delete(key);
      });
    }
    return pending;
  } catch (error) {
    return Promise.reject(error);
  }
};

/**
 * @param {unknown} file A `File`, a Session resource view, or any value (rejected unless it has
 *   `arrayBuffer()` or `text()`): each member is tested before it is used.
 * @returns {Promise<string>}
 */
export const readFileText = async (file) => {
  const lazyText = isSessionResourceFileView(file)
    ? readSessionResourceText(file)
    : null;
  if (lazyText) return lazyText;
  // A view: the reads below test each member they use.
  const view = /** @type {{ readonly arrayBuffer?: unknown, readonly text?: unknown } | null | undefined} */ (file);
  if (typeof view?.arrayBuffer === 'function') {
    return bytesToText(await readFileBytes(file));
  }
  // A function member is a File's or Blob's (Session JSON holds none), whose `text()` resolves to a string.
  if (typeof view?.text === 'function') return /** @type {{ text: () => Promise<string> }} */ (file).text();
  throw new TypeError('A File-like object with arrayBuffer() or text() is required.');
};

/**
 * @param {Uint8Array | ArrayBuffer} bytes
 * @returns {Promise<ResourceContentToken>}
 */
export const resourceContentToken = async (bytes) => /** @type {ResourceContentToken} */ (await sha256Hex(bytes));

/**
 * The content token of a File or Session resource view, hashed once per object.
 * @param {FileLike} file
 * @returns {Promise<ResourceContentToken>}
 */
export const readFileContentToken = (file) => {
  let pending = fileContentTokenCache.get(file);
  if (!pending) {
    pending = readFileBytes(file).then(resourceContentToken);
    fileContentTokenCache.set(file, pending);
    pending.catch(() => {
      if (fileContentTokenCache.get(file) === pending) fileContentTokenCache.delete(file);
    });
  }
  return pending;
};

/** @param {unknown} file */
export const getSessionResourceSource = (file) => sessionResourceSource(file);

/** @param {FileLike | null | undefined} file */
export const releaseFileContentCache = (file) => {
  const key = asCacheKey(file);
  if (!key) return false;
  const deleted = fileByteCache.delete(key);
  return releaseSessionResourceContent(file) || deleted;
};

/** @param {FileLike} file */
export const takeFileBytesForTransfer = async (file) => {
  const bytes = await readFileBytes(file);
  const buffer = (
    bytes.byteOffset === 0 && bytes.byteLength === bytes.buffer.byteLength
      ? bytes.buffer
      : bytes.slice().buffer
  );
  releaseFileContentCache(file);
  return buffer;
};

/** @param {FileLike} file */
export const cloneFileBytesForTransfer = async (file) => {
  const bytes = await readFileBytes(file);
  return bytes.slice().buffer;
};
