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
  expandGzipBytes,
  sha256Digest,
  sha256Hex,
  textToBase64,
  textToBytes
} from './byte-utils.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { MAX_EXPANDED_INPUT_BYTES } from '../contracts/schema.generated.js';
import { testByteCap } from './runtime-test-hooks.js';

export {
  base64ToBytes,
  bytesToBase64,
  bytesToText,
  sha256Hex,
  textToBase64,
  textToBytes
};

/** @import { FileLike, RawSessionResourceDescriptor } from './session-resource-backing.js' */

/**
 * A resource's content token: the Worker's cache key for its bytes, so equal
 * bytes share one token whatever resource IDs name them (`contentTokenOfBytes`).
 * Only the page and the Worker compare tokens: Python checks their form, and no
 * Session holds one.
 * @typedef {string & { readonly __brand: 'ResourceContentToken' }} ResourceContentToken
 */
/**
 * What the page knows of a File's content from its first read: its token and
 * its byte length (gzip expanded). A File's bytes never change, so these
 * outlive the release of its bytes.
 * @typedef {{ token: ResourceContentToken, byteLength: number }} FileContentFacts
 */

/** @type {WeakMap<object, Promise<Uint8Array>>} */
const fileByteCache = new WeakMap();
/** @type {WeakMap<object, Promise<FileContentFacts>>} */
const fileContentFacts = new WeakMap();
// The Files whose cached bytes a gzip file expanded to. They stay until a
// transfer takes them: expanding again would be slower than keeping them.
/** @type {WeakSet<object>} */
const expandedFiles = new WeakSet();
/**
 * The descriptor Generate made of each File (`fileBackedDescriptor`), keyed by that File. A
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

// The DOMException names a browser gives when it refuses to read a File (File
// API, "file read errors"): the file was moved, deleted, or changed on disk
// after it was picked, or the browser deems it unsafe to read.
const REFUSED_FILE_READS = new Set(['NotFoundError', 'NotReadableError', 'SecurityError']);
/** @param {unknown} error */
const refusedReadError = (error) => (
  error instanceof DOMException && REFUSED_FILE_READS.has(error.name)
    ? diagnosticError('INPUT_UNREADABLE', {}, { stage: 'read' })
    : error
);

const MIB = 1024 * 1024;
const TOKEN_PIECE_BYTES = 16 * MIB;

/**
 * The content token of `bytes`, one scheme for every size: the SHA-256 hex of
 * the SHA-256 digests of their 16 MiB pieces in order (one piece when empty).
 * One piece is hashed at a time, so hashing copies at most one piece.
 * @param {Uint8Array} bytes
 * @returns {Promise<ResourceContentToken>}
 */
export const contentTokenOfBytes = async (bytes) => {
  const pieceCount = Math.max(1, Math.ceil(bytes.byteLength / TOKEN_PIECE_BYTES));
  const digests = new Uint8Array(pieceCount * 32);
  for (let index = 0; index < pieceCount; index += 1) {
    const start = index * TOKEN_PIECE_BYTES;
    digests.set(new Uint8Array(await sha256Digest(bytes.subarray(start, start + TOKEN_PIECE_BYTES))), index * 32);
  }
  return /** @type {ResourceContentToken} */ (await sha256Hex(digests));
};

/**
 * The facts of `key`, started from `bytes` unless an earlier read started them.
 * @param {object} key @param {Promise<Uint8Array>} bytes
 * @returns {Promise<FileContentFacts>}
 */
const startContentFacts = (key, bytes) => {
  let facts = fileContentFacts.get(key);
  if (!facts) {
    facts = bytes.then(async (read) => ({ token: await contentTokenOfBytes(read), byteLength: read.byteLength }));
    fileContentFacts.set(key, facts);
    facts.catch(() => {
      if (fileContentFacts.get(key) === facts) fileContentFacts.delete(key);
    });
  }
  return facts;
};

// A File's bytes are kept while a hold covers its read (a Generate, a Save, a
// once-per-File read): the reads it covers share them. When the last hold ends,
// or a read settles with no hold, the page drops every settled plain read; the
// next read reads the File again.
let readHolds = 0;
/** @type {Map<object, Promise<Uint8Array>>} */
const sharedReads = new Map();
/** @type {WeakSet<Promise<Uint8Array>>} */
const settledReads = new WeakSet();
const dropSharedReads = () => {
  if (readHolds) return;
  for (const [key, pending] of sharedReads) {
    if (!settledReads.has(pending)) continue;
    sharedReads.delete(key);
    if (fileByteCache.get(key) === pending && !expandedFiles.has(key)) fileByteCache.delete(key);
  }
};

/**
 * Keeps the bytes of the Files read until the returned release is called.
 * @returns {() => void} Releases this hold; a second call does nothing.
 */
export const holdReadBytes = () => {
  readHolds += 1;
  let held = true;
  return () => {
    if (!held) return;
    held = false;
    readHolds -= 1;
    dropSharedReads();
  };
};
/**
 * The plain bytes of a File's bytes: a gzip file expands to at most Python's
 * `MAX_EXPANDED_INPUT_BYTES`, and past it rejects with `INPUT_UNREADABLE`
 * (`GZIP_EXPANDED_SIZE`), naming the file and the cap in MiB.
 * @param {string} name The File's name.
 * @param {Uint8Array<ArrayBuffer>} bytes
 */
const expandInputBytes = (name, bytes) => {
  const maxBytes = testByteCap('maxExpandedInputBytes') ?? MAX_EXPANDED_INPUT_BYTES;
  return expandGzipBytes(bytes, maxBytes, () => diagnosticError('INPUT_UNREADABLE',
    { reason: 'GZIP_EXPANDED_SIZE', inputName: name, limit: Math.floor(maxBytes / MIB) }, { stage: 'read' }));
};

/**
 * Reads `file` (cached under `key`), starts its content facts, and shares the
 * bytes with the reads the current holds cover (`holdReadBytes`).
 * @param {object} key
 * @param {{ arrayBuffer: () => Promise<ArrayBuffer>, readonly name?: unknown }} file
 * @returns {Promise<Uint8Array>}
 */
const startFileRead = (key, file) => {
  const pending = Promise.resolve()
    .then(() => file.arrayBuffer())
    .then(asBytes, (error) => { throw refusedReadError(error); })
    .then(async (bytes) => {
      const plain = await expandInputBytes(String(file.name ?? ''), bytes);
      if (plain !== bytes) expandedFiles.add(key);
      return plain;
    });
  fileByteCache.set(key, pending);
  sharedReads.set(key, pending);
  startContentFacts(key, pending);
  const settle = () => {
    settledReads.add(pending);
    dropSharedReads();
  };
  pending.then(settle, () => {
    if (fileByteCache.get(key) === pending) fileByteCache.delete(key);
    settle();
  });
  return pending;
};

/**
 * The bytes of an uploaded File, the one place they are read: a gzip file
 * (NCBI's `.gbff.gz`, `.fna.gz`, `.gff.gz`) reads as the plain file it holds,
 * told by its first two bytes, so every browser reader, the Worker and the
 * Session see the plain bytes. A read the browser refuses rejects with
 * `INPUT_UNREADABLE`, as a gzip file over the expanded cap does.
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
    // The `typeof` check above holds its `arrayBuffer()`.
    return fileByteCache.get(key)
      ?? startFileRead(key, /** @type {{ arrayBuffer: () => Promise<ArrayBuffer>, readonly name?: unknown }} */ (file));
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
 * The content facts of a File or Session resource view, computed once per
 * object from its first read; known facts read nothing.
 * @param {FileLike} file
 * @returns {Promise<FileContentFacts>}
 */
export const readFileContentFacts = (file) => {
  const key = asCacheKey(file);
  if (!key) return Promise.reject(new TypeError('A File-like object is required.'));
  return fileContentFacts.get(key) ?? startContentFacts(key, readFileBytes(file));
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
 * @returns {RawSessionResourceDescriptor | null} The descriptor Generate made of `file`.
 */
export const fileDescriptorOf = (file) => {
  const key = asCacheKey(file);
  return (key && fileDescriptors.get(key)) || null;
};

/** @param {FileLike | null | undefined} file */
export const releaseFileContentCache = (file) => {
  const key = asCacheKey(file);
  if (!key) return false;
  sharedReads.delete(key);
  const deleted = fileByteCache.delete(key);
  return releaseSessionResourceContent(file) || deleted;
};

// The buffers a take handed out whole. Each is its caller's to transfer: a take
// that read the same cached bytes (it started before an earlier take released
// them, as two helper requests staging one input at once do, OV-2490) reads
// the file again instead of sharing a buffer that the first transfer detaches.
// A read again that returns a taken buffer too (a FileLike that keeps one
// buffer) gives a copy, so the retake ends.
/** @type {WeakSet<ArrayBufferLike>} */
const takenBuffers = new WeakSet();

/** @param {FileLike} file */
const readAndRelease = async (file) => {
  const bytes = await readFileBytes(file);
  releaseFileContentCache(file);
  return bytes;
};

/**
 * @param {FileLike} file
 * @returns {Promise<ArrayBuffer>} A buffer no other take holds.
 */
export const takeFileBytesForTransfer = async (file) => {
  let bytes = await readAndRelease(file);
  if (takenBuffers.has(bytes.buffer)) bytes = await readAndRelease(file);
  if (takenBuffers.has(bytes.buffer) || bytes.byteOffset !== 0 || bytes.byteLength !== bytes.buffer.byteLength) {
    return bytes.slice().buffer;
  }
  takenBuffers.add(bytes.buffer);
  return /** @type {ArrayBuffer} */ (bytes.buffer);
};

/** @param {FileLike} file */
export const cloneFileBytesForTransfer = async (file) => {
  const bytes = await readFileBytes(file);
  return bytes.slice().buffer;
};
