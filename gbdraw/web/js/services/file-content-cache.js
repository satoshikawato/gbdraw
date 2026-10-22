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

/** @import { FileLike, RawSessionResourceDescriptor } from './session-resource-backing.js' */

/** @type {WeakMap<object, Promise<Uint8Array>>} */
const fileByteCache = new WeakMap();
/**
 * The descriptor each File was encoded or loaded as, keyed by that File. A
 * File never changes its bytes, so the File object is the key.
 * @type {WeakMap<object, RawSessionResourceDescriptor>}
 */
const fileDescriptors = new WeakMap();

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

/** @param {unknown} file */
export const getSessionResourceSource = (file) => sessionResourceSource(file);

/**
 * @template {RawSessionResourceDescriptor} D
 * @param {object} file A `File`, or another object whose bytes were read.
 * @param {D} descriptor
 * @returns {D}
 */
export const rememberFileDescriptor = (file, descriptor) => {
  fileDescriptors.set(file, descriptor);
  return descriptor;
};

/**
 * @param {unknown} file
 * @returns {RawSessionResourceDescriptor | null} The descriptor `file` was encoded or loaded as.
 */
export const fileDescriptorOf = (file) => {
  const key = asCacheKey(file);
  return (key && fileDescriptors.get(key)) || null;
};

/**
 * The base64 text a File was encoded or loaded as, when its byte size matches.
 * @typedef {Readonly<{ encoding: 'base64', data: string, size: number }>} EncodedFilePayload
 */
/**
 * @param {FileLike} file
 * @returns {EncodedFilePayload | null}
 */
export const encodedFilePayloadOf = (file) => {
  const descriptor = fileDescriptorOf(file);
  const size = descriptor?.size;
  if (
    !descriptor || descriptor.encoding !== 'base64' || typeof descriptor.data !== 'string'
    || typeof size !== 'number' || size !== file.size
  ) return null;
  return Object.freeze({ encoding: 'base64', data: descriptor.data, size });
};

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
