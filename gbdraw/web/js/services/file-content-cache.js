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
  expandGzipPieces,
  expandsAsGzip,
  joinBytes,
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
/** @import { InputCopy, InputCopyStore, InputCopyWriter } from './input-copies.js' */

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
// The Files whose cached bytes a gzip file expanded to while no copy holds
// them. They stay until a transfer takes them: expanding again would be slower
// than keeping them.
/** @type {WeakSet<object>} */
const expandedFiles = new WeakSet();
// The temporary copies of the chosen inputs (PD-OI-121, `services/input-copies.js`):
// the store the composition root gives (`useInputCopyStore`), the Files the
// input store holds (`keepFileCopy`), whose first read writes a copy, the copy
// of each such File, and the Files the input store let go, whose copies are
// given back once no reading runs (`releaseFileCopy`).
/** @type {InputCopyStore | null} */
let copyStore = null;
/** @type {WeakSet<object>} */
const keptFiles = new WeakSet();
/** @type {WeakMap<object, Promise<InputCopy | null>>} */
const fileCopies = new WeakMap();
/** @type {Set<object>} */
const releasedFiles = new Set();
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
 * @typedef {object} PieceSink
 * @property {(piece: Uint8Array, owned: boolean) => Promise<void>} add Hashes one piece and starts its write.
 * @property {() => Promise<{ facts: FileContentFacts, copy: InputCopy | null }>} finish
 * @property {() => void} fail Deletes the copy.
 */
/**
 * Hashes the pieces of one content in order and, given a writer, writes each
 * piece to the content's copy as it is hashed: a piece is hashed while the one
 * before it is written, so at most one piece waits for the Worker. A copy that
 * fails to write is deleted, and the content keeps no copy.
 * @param {InputCopyWriter | null} writer
 * @returns {PieceSink}
 */
const createPieceSink = (writer) => {
  /** @type {Uint8Array[]} */
  const digests = [];
  let byteLength = 0;
  let copying = writer;
  /** @type {Promise<void>} */
  let writing = Promise.resolve();
  const stopCopying = () => {
    copying?.abort();
    copying = null;
  };
  return {
    // `owned`: the piece's buffer is the sink's to transfer.
    add: async (piece, owned) => {
      // Web Crypto copies the piece when it is called, so the piece may then be transferred.
      const digest = sha256Digest(piece);
      byteLength += piece.byteLength;
      if (copying) {
        const whole = owned && piece.byteOffset === 0 && piece.byteLength === piece.buffer.byteLength;
        const buffer = /** @type {ArrayBuffer} */ (whole ? piece.buffer : piece.slice().buffer);
        await writing;
        if (copying) writing = copying.write(buffer).catch(stopCopying);
      }
      digests.push(new Uint8Array(await digest));
    },
    finish: async () => {
      await writing;
      const token = /** @type {ResourceContentToken} */ (await sha256Hex(joinBytes(digests)));
      const copy = copying ? await copying.finish(token).catch(() => null) : null;
      return { facts: { token, byteLength }, copy };
    },
    fail: stopCopying
  };
};

/**
 * Hands `bytes` to `sink` in 16 MiB pieces (one empty piece when empty).
 * @param {Uint8Array} bytes @param {PieceSink} sink
 */
const addPlainPieces = async (bytes, sink) => {
  for (let start = 0; start === 0 || start < bytes.byteLength; start += TOKEN_PIECE_BYTES) {
    await sink.add(bytes.subarray(start, start + TOKEN_PIECE_BYTES), false);
  }
};

/**
 * The content token of `bytes`, one scheme for every size: the SHA-256 hex of
 * the SHA-256 digests of their 16 MiB pieces in order (one piece when empty).
 * One piece is hashed at a time, so hashing copies at most one piece.
 * @param {Uint8Array} bytes
 * @returns {Promise<ResourceContentToken>}
 */
export const contentTokenOfBytes = async (bytes) => {
  const sink = createPieceSink(null);
  await addPlainPieces(bytes, sink);
  return (await sink.finish()).facts.token;
};

/**
 * Keeps `facts` as the facts of `key` until they fail.
 * @param {object} key @param {Promise<FileContentFacts>} facts
 */
const rememberFacts = (key, facts) => {
  fileContentFacts.set(key, facts);
  facts.catch(() => {
    if (fileContentFacts.get(key) === facts) fileContentFacts.delete(key);
  });
  return facts;
};

// A File's bytes are kept while a hold covers its read (a Generate, a Save, a
// once-per-File read): the reads it covers share them. When the last hold ends,
// or a read settles with no hold, the page drops every settled plain read; the
// next read reads the File's copy, or the File again. A first read settles once
// its copy is written, so a read that starts meanwhile shares its bytes.
let readHolds = 0;
/** @type {Map<object, Promise<Uint8Array>>} */
const sharedReads = new Map();
/** @type {WeakSet<Promise<Uint8Array>>} */
const settledReads = new WeakSet();
/** @param {object} key */
const readRuns = (key) => {
  const pending = sharedReads.get(key);
  return Boolean(pending && !settledReads.has(pending));
};
const dropSharedReads = () => {
  if (readHolds) return;
  for (const [key, pending] of sharedReads) {
    if (!settledReads.has(pending)) continue;
    sharedReads.delete(key);
    if (fileByteCache.get(key) === pending && !expandedFiles.has(key)) fileByteCache.delete(key);
  }
  for (const key of releasedFiles) {
    if (readRuns(key)) continue;
    releasedFiles.delete(key);
    fileCopies.get(key)?.then((copy) => copy?.release());
    fileCopies.delete(key);
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

const expandedCap = () => testByteCap('maxExpandedInputBytes') ?? MAX_EXPANDED_INPUT_BYTES;
/**
 * The refusal of a gzip file that expands past Python's `MAX_EXPANDED_INPUT_BYTES`:
 * `INPUT_UNREADABLE` (`GZIP_EXPANDED_SIZE`), naming the file and the cap in MiB.
 * @param {string} name @param {number} maxBytes
 */
const expandedSizeError = (name, maxBytes) => diagnosticError('INPUT_UNREADABLE',
  { reason: 'GZIP_EXPANDED_SIZE', inputName: name, limit: Math.floor(maxBytes / MIB) }, { stage: 'read' });

/**
 * The plain bytes of a File's bytes: a gzip file expands, up to the cap.
 * @param {string} name The File's name.
 * @param {Uint8Array<ArrayBuffer>} bytes
 * @returns {Promise<Uint8Array>}
 */
const expandInputBytes = async (name, bytes) => {
  if (!expandsAsGzip(bytes)) return bytes;
  /** @type {Uint8Array[]} */
  const pieces = [];
  const maxBytes = expandedCap();
  await expandGzipPieces(bytes, maxBytes, () => expandedSizeError(name, maxBytes), TOKEN_PIECE_BYTES, async (piece) => {
    pieces.push(piece);
  });
  return joinBytes(pieces);
};

/** @typedef {{ arrayBuffer: () => Promise<ArrayBuffer>, readonly name?: unknown }} ReadableFile */

/**
 * The File's own bytes, as stored: a read the browser refuses rejects with `INPUT_UNREADABLE`.
 * @param {ReadableFile} file
 */
const readSource = (file) => Promise.resolve()
  .then(() => file.arrayBuffer())
  .then(asBytes, (error) => { throw refusedReadError(error); });

/**
 * A later read of a File that has no copy: its bytes are checked against the
 * facts of its first read, and a File whose bytes changed on disk since (a
 * browser may return the new bytes instead of refusing) rejects with
 * `INPUT_UNREADABLE`, as a refused read does (PD-OI-121).
 * @param {ReadableFile} file @param {Promise<FileContentFacts>} knownFacts
 */
const readAgain = async (file, knownFacts) => {
  const bytes = await expandInputBytes(String(file.name ?? ''), await readSource(file));
  const facts = await knownFacts;
  if (bytes.byteLength !== facts.byteLength || await contentTokenOfBytes(bytes) !== facts.token) {
    throw diagnosticError('INPUT_UNREADABLE', {}, { stage: 'read' });
  }
  return bytes;
};

/**
 * A later read through a File's copy, or of the File itself when it has no copy or the copy is gone.
 * @param {ReadableFile} file @param {Promise<InputCopy | null>} copy @param {Promise<FileContentFacts>} facts
 */
const readThroughCopy = async (file, copy, facts) => {
  const held = await copy;
  if (held) {
    try {
      return await held.read();
    } catch {
      // The copy is gone (the site's storage was cleared): the File is read itself.
    }
  }
  return readAgain(file, facts);
};

/**
 * The first read of a File: the one read of its own bytes. It hashes them in
 * 16 MiB pieces (a gzip file's as they expand) for the File's facts and, for a
 * File the input store holds, writes each piece to the File's copy. Readers get
 * a plain File's bytes as soon as they are read; a gzip file's expanded bytes
 * come from its copy, so no whole expanded buffer is held, or from its pieces
 * when it has no copy.
 * @param {object} key @param {ReadableFile} file
 * @returns {{ bytes: Promise<Uint8Array>, done: Promise<unknown> }}
 */
const readFirst = (key, file) => {
  const source = readSource(file);
  const written = source.then(async (raw) => {
    const writer = keptFiles.has(key) && copyStore ? await copyStore.openWriter().catch(() => null) : null;
    const sink = createPieceSink(writer);
    /** @type {Uint8Array[] | null} */
    const pieces = writer ? null : [];
    try {
      if (!expandsAsGzip(raw)) await addPlainPieces(raw, sink);
      else {
        const maxBytes = expandedCap();
        const name = String(file.name ?? '');
        await expandGzipPieces(raw, maxBytes, () => expandedSizeError(name, maxBytes), TOKEN_PIECE_BYTES, async (piece) => {
          pieces?.push(piece);
          await sink.add(piece, true);
        });
      }
    } catch (error) {
      sink.fail();
      throw error;
    }
    return { ...await sink.finish(), pieces };
  });
  const facts = rememberFacts(key, written.then((result) => result.facts));
  const copy = written.then((result) => result.copy, () => null);
  fileCopies.set(key, copy);
  written.catch(() => {
    if (fileCopies.get(key) === copy) fileCopies.delete(key);
  });
  const bytes = source.then((raw) => (expandsAsGzip(raw) ? written.then(({ pieces }) => {
    if (!pieces) return readThroughCopy(file, copy, facts);
    expandedFiles.add(key);
    return joinBytes(pieces);
  }) : raw));
  return { bytes, done: Promise.allSettled([bytes, written]) };
};

/**
 * Reads `file` (cached under `key`): its first read, or a later read through
 * its copy, and shares the bytes with the reads the current holds cover
 * (`holdReadBytes`).
 * @param {object} key
 * @param {ReadableFile} file
 * @returns {Promise<Uint8Array>}
 */
const startFileRead = (key, file) => {
  const facts = fileContentFacts.get(key);
  const read = facts
    ? { bytes: readThroughCopy(file, fileCopies.get(key) ?? Promise.resolve(null), facts) }
    : readFirst(key, file);
  const pending = read.bytes;
  fileByteCache.set(key, pending);
  sharedReads.set(key, pending);
  pending.catch(() => {
    if (fileByteCache.get(key) === pending) fileByteCache.delete(key);
  });
  const settle = () => {
    settledReads.add(pending);
    dropSharedReads();
  };
  ('done' in read ? read.done : pending).then(settle, settle);
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
    return fileByteCache.get(key) ?? startFileRead(key, /** @type {ReadableFile} */ (file));
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
  const known = fileContentFacts.get(key);
  if (known) return known;
  // A File's first read keeps its facts; a Session resource view's are hashed from its bytes.
  const bytes = readFileBytes(file);
  return fileContentFacts.get(key) ?? rememberFacts(key, bytes.then(async (read) => ({
    token: await contentTokenOfBytes(read), byteLength: read.byteLength
  })));
};

/**
 * The store the temporary copies are written to; null writes none, and every File is read itself.
 * @param {InputCopyStore | null} store
 */
export const useInputCopyStore = (store) => {
  copyStore = store;
};

/**
 * The input store now holds `file`: the File is read and copied now (PD-OI-121),
 * unless a read of it has started, and every later reading reads the copy.
 * @param {unknown} file
 */
export const keepFileCopy = (file) => {
  const key = asCacheKey(file);
  if (!key || !copyStore || isSessionResourceFileView(file)) return;
  if (typeof /** @type {{ readonly arrayBuffer?: unknown }} */ (key).arrayBuffer !== 'function') return;
  releasedFiles.delete(key);
  if (keptFiles.has(key)) return;
  keptFiles.add(key);
  if (!fileContentFacts.has(key)) readFileContentFacts(/** @type {FileLike} */ (file)).catch(() => {});
};

/**
 * The input store let go of `file`: its copy is deleted once no hold and no read of it runs.
 * @param {unknown} file
 */
export const releaseFileCopy = (file) => {
  const key = asCacheKey(file);
  if (!key || !keptFiles.has(key)) return;
  keptFiles.delete(key);
  releasedFiles.add(key);
  dropSharedReads();
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
  // A first read hashes and copies the bytes it read until its facts are known; a transfer detaches them after that.
  const key = asCacheKey(file);
  await (key && fileContentFacts.get(key))?.catch(() => {});
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
