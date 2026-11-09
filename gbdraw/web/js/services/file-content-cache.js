// @ts-check
import {
  checkInputExpansion,
  expandInputBytes,
  isInputExpansionRefusal,
  isSessionResourceFileView,
  releaseSessionResourceContent,
  readSessionResourceBytes,
  readSessionResourceStoredBytes,
  readSessionResourceText,
  sessionResourceSource,
  sessionResourceStoredBytes
} from './session-resource-backing.js';
import {
  asBytes,
  base64ToBytes,
  bytesToBase64,
  bytesToText,
  sha256Digest,
  sha256Hex,
  textToBase64,
  textToBytes
} from './byte-utils.js';
import { diagnosticError, normalizeUserFacingError } from '../utils/error-normalization.js';

export {
  base64ToBytes,
  bytesToBase64,
  bytesToText,
  sha256Hex,
  textToBase64,
  textToBytes
};

/** @import { FileLike, RawSessionResourceDescriptor, SessionResourceFileView, StoredResourceBytes } from './session-resource-backing.js' */
/** @import { InputCopy, InputCopyStore, InputCopyWriter } from './input-copies.js' */
/** @import { WorkerInputFile } from './diagram-worker-protocol.js' */

/**
 * A resource's content token: the Worker's cache key for its bytes, so equal
 * bytes share one token whatever resource IDs name them (`contentTokenOfBytes`).
 * Only the page and the Worker compare tokens: Python checks their form, and no
 * Session holds one.
 * @typedef {string & { readonly __brand: 'ResourceContentToken' }} ResourceContentToken
 */
/**
 * What the page knows of a content from its first read: the token and byte
 * length of its bytes as chosen (a File) or as stored (a Session resource
 * view), a gzip input's gzip bytes. Its copy holds them, Generate sends them
 * (Python expands a gzip input), and a Session file stores them (receipt
 * session.file-container). A File's bytes never change, so these outlive the
 * release of its bytes.
 * @typedef {{ token: ResourceContentToken, byteLength: number }} ContentFacts
 */

// The bytes as chosen of each File read, kept while a hold covers them
// (`holdReadBytes`); a gzip file's expanded bytes are never kept.
/** @type {WeakMap<object, Promise<Uint8Array>>} */
const fileByteCache = new WeakMap();
/** @type {WeakMap<object, Promise<ContentFacts>>} */
const fileContentFacts = new WeakMap();
// The temporary copies of the chosen inputs (PD-OI-121, `services/input-copies.js`):
// the store `state.js` gives (`useInputCopyStore`), the Files the input store
// holds (`keepFileCopy`), whose first read writes a copy, the other uses that
// keep a copy (`holdFileCopy`), the copy of each such File, and the Files no
// longer kept or used, whose copies are given back once no reading runs.
/** @type {InputCopyStore | null} */
let copyStore = null;
/** @type {WeakSet<object>} */
const keptFiles = new WeakSet();
/** @type {WeakMap<object, Promise<InputCopy | null>>} */
const fileCopies = new WeakMap();
/** @type {WeakMap<object, number>} */
const copyUses = new WeakMap();
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
 * The content token of `bytes`, one scheme for every size: the SHA-256 hex of
 * the SHA-256 digests of their 16 MiB pieces in order (one piece when empty).
 * One piece is hashed at a time, so hashing copies at most one piece.
 * @param {Uint8Array} bytes
 * @returns {Promise<ResourceContentToken>}
 */
export const contentTokenOfBytes = async (bytes) => {
  const pieceCount = Math.max(1, Math.ceil(bytes.byteLength / TOKEN_PIECE_BYTES));
  /** @type {ArrayBuffer[]} */
  const digests = [];
  for (let index = 0; index < pieceCount; index += 1) {
    const start = index * TOKEN_PIECE_BYTES;
    digests.push(await sha256Digest(bytes.subarray(start, start + TOKEN_PIECE_BYTES)));
  }
  return tokenOfPieceDigests(digests);
};

/**
 * The content token from the SHA-256 digests of the 16 MiB pieces, in order (`contentTokenOfBytes`).
 * @param {readonly ArrayBuffer[]} digests
 * @returns {Promise<ResourceContentToken>}
 */
const tokenOfPieceDigests = async (digests) => {
  const joined = new Uint8Array(digests.length * 32);
  digests.forEach((digest, index) => joined.set(new Uint8Array(digest), index * 32));
  return /** @type {ResourceContentToken} */ (await sha256Hex(joined));
};

/**
 * The facts of `bytes`, their length taken before they are hashed.
 * @param {Uint8Array} bytes
 * @returns {Promise<ContentFacts>}
 */
const factsOfBytes = async (bytes) => {
  const byteLength = bytes.byteLength;
  return { token: await contentTokenOfBytes(bytes), byteLength };
};

/**
 * Keeps `facts` as the facts of `key` until they fail.
 * @param {object} key @param {Promise<ContentFacts>} facts
 */
const rememberFacts = (key, facts) => {
  fileContentFacts.set(key, facts);
  facts.catch(() => {
    if (fileContentFacts.get(key) === facts) fileContentFacts.delete(key);
  });
  return facts;
};

/**
 * Cuts bytes, as they come, into the 16 MiB pieces of the content token and
 * hashes each (`contentTokenOfBytes`), the one owner of a streamed token. Each
 * piece goes to `write` too when given, which may transfer it and give a
 * buffer back for the next piece, so the pieces go through one buffer.
 * @param {((piece: ArrayBuffer) => Promise<ArrayBuffer>) | null} [write]
 */
const pieceHasher = (write = null) => {
  /** @type {ArrayBuffer[]} */
  const digests = [];
  let byteLength = 0;
  let buffer = new ArrayBuffer(TOKEN_PIECE_BYTES);
  let filled = 0;
  // Web Crypto copies a piece when the digest starts, so the write may transfer it at once.
  const flush = async () => {
    const piece = filled === buffer.byteLength ? buffer : buffer.slice(0, filled);
    const [digest, given] = await Promise.all([sha256Digest(new Uint8Array(piece)), write ? write(piece) : null]);
    digests.push(digest);
    byteLength += filled;
    if (piece === buffer && given) buffer = given.byteLength === TOKEN_PIECE_BYTES ? given : new ArrayBuffer(TOKEN_PIECE_BYTES);
    filled = 0;
  };
  return {
    /** @param {Uint8Array} chunk */
    async add(chunk) {
      for (let rest = chunk; rest.byteLength > 0;) {
        const taken = Math.min(rest.byteLength, buffer.byteLength - filled);
        new Uint8Array(buffer, filled, taken).set(rest.subarray(0, taken));
        filled += taken;
        rest = rest.subarray(taken);
        if (filled === buffer.byteLength) await flush();
      }
    },
    // The last piece, or the one empty piece of empty content (`contentTokenOfBytes`).
    /** @returns {Promise<ContentFacts>} */
    async finish() {
      if (filled > 0 || digests.length === 0) await flush();
      return { token: await tokenOfPieceDigests(digests), byteLength };
    }
  };
};

/**
 * The facts of a stream's bytes, read once and one piece at a time: the
 * content is never whole in memory.
 * @param {ReadableStream<Uint8Array>} stream
 * @returns {Promise<ContentFacts>}
 */
export const contentFactsOfStream = async (stream) => {
  const hasher = pieceHasher();
  const reader = stream.getReader();
  try {
    for (let step = await reader.read(); !step.done; step = await reader.read()) await hasher.add(step.value);
  } catch (error) {
    reader.cancel().catch(() => {});
    throw error;
  }
  return hasher.finish();
};

// A File's bytes are kept while a hold covers its read (a Generate, a Save, a
// once-per-File read): the reads it covers share them. When the last hold ends,
// or a read settles with no hold, the page drops every settled read; the
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
    if (fileByteCache.get(key) === pending) fileByteCache.delete(key);
  }
  for (const key of releasedFiles) {
    // A File kept or used again keeps its copy; one that a reading still reads gives it back after.
    if (keptFiles.has(key) || copyUses.get(key)) releasedFiles.delete(key);
    if (readRuns(key) || !releasedFiles.has(key)) continue;
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

/**
 * A File as the cache reads it; a real File also streams and has a size (a test's stand-in may not).
 * @typedef {{
 *   arrayBuffer: () => Promise<ArrayBuffer>, readonly name?: unknown,
 *   readonly stream?: () => ReadableStream<Uint8Array>, readonly size?: unknown
 * }} ReadableFile
 */

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
 * @param {ReadableFile} file @param {Promise<ContentFacts>} knownFacts
 */
const readAgain = async (file, knownFacts) => {
  const bytes = await readSource(file);
  const facts = await knownFacts;
  if (bytes.byteLength !== facts.byteLength || await contentTokenOfBytes(bytes) !== facts.token) {
    throw diagnosticError('INPUT_UNREADABLE', {}, { stage: 'read' });
  }
  return bytes;
};

/**
 * A later read of a File's bytes as chosen: through its copy, or of the File
 * itself when it has no copy or the copy is gone.
 * @param {ReadableFile} file @param {Promise<InputCopy | null>} copy @param {Promise<ContentFacts>} facts
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
 * Writes a File's bytes as chosen (`raw`) to a new copy in 16 MiB pieces, one
 * at a time, while their facts are hashed, and closes it under their token,
 * so a copy is shared only by Files and Session entries of the same bytes;
 * null when no copy could be written, or when it holds fewer bytes than the
 * facts count. Every piece goes through one buffer, which the Worker gives
 * back, so the copy holds one piece beside `raw`.
 * @param {InputCopyStore} store @param {Uint8Array} raw @param {Promise<ContentFacts>} facts
 * @returns {Promise<InputCopy | null>}
 */
const writeCopy = async (store, raw, facts) => {
  const writer = await store.openWriter().catch(() => null);
  if (!writer) return null;
  try {
    /** @type {ArrayBuffer | null} */
    let spare = null;
    let written = 0;
    for (let start = 0; start === 0 || start < raw.byteLength; start += TOKEN_PIECE_BYTES) {
      const piece = raw.subarray(start, start + TOKEN_PIECE_BYTES);
      const buffer = spare?.byteLength === piece.byteLength ? spare : new ArrayBuffer(piece.byteLength);
      new Uint8Array(buffer).set(piece);
      spare = await writer.write(buffer);
      written += piece.byteLength;
    }
    const { token, byteLength } = await facts;
    if (written !== byteLength) throw new Error(`Copied ${written} of ${byteLength} bytes.`);
    return await writer.finish(token);
  } catch {
    writer.abort();
    return null;
  }
};

/**
 * A copy of a stream's bytes and their facts; or no copy, because the browser
 * lets none be written (`unavailable`), or a write failed or fell short, as when
 * the browser's storage is full (`write-failed`). A failed copy is deleted.
 * @typedef {{ copy: InputCopy, facts: ContentFacts }
 *   | { copy: null, reason: 'unavailable' | 'write-failed' }} StreamCopy
 */
/** A write of a stream's copy that failed, told apart from a failure of the stream. */
const COPY_WRITE_FAILED = Symbol('copy write failed');

/**
 * Copies the bytes of `stream` (an input that a Session Load reads from its
 * file, PD-OI-121) into a new copy in 16 MiB pieces, one at a time, hashing each
 * piece for the content token while it is written: the content is never whole
 * in memory and is read once. A copy that holds the token already is shared.
 * A stream that fails rejects with its error, after the copy is deleted; with
 * no copy store, the stream is cancelled unread.
 * @param {ReadableStream<Uint8Array>} stream
 * @returns {Promise<StreamCopy>}
 */
export const copyInputStream = async (stream) => {
  const writer = await openCopyWriter();
  if (!writer) {
    await stream.cancel().catch(() => {});
    return { copy: null, reason: 'unavailable' };
  }
  return writeStreamCopy(writer, stream);
};

/** Whether the page has a copy store (`useInputCopyStore`); its writes may still fail. */
export const inputCopiesAvailable = () => copyStore !== null;

/** A writer of a new copy, or null when the browser lets none be written. */
const openCopyWriter = async () => (copyStore ? copyStore.openWriter().catch(() => null) : null);

/**
 * Writes `stream` through `writer` as `copyInputStream` describes.
 * @param {InputCopyWriter} writer @param {ReadableStream<Uint8Array>} stream
 * @returns {Promise<StreamCopy>}
 */
const writeStreamCopy = async (writer, stream) => {
  const reader = stream.getReader();
  const hasher = pieceHasher((piece) => writer.write(piece).catch(() => { throw COPY_WRITE_FAILED; }));
  try {
    for (let step = await reader.read(); !step.done; step = await reader.read()) await hasher.add(step.value);
    const facts = await hasher.finish();
    const copy = await writer.finish(facts.token).catch(() => { throw COPY_WRITE_FAILED; });
    return { copy, facts };
  } catch (error) {
    writer.abort();
    reader.cancel().catch(() => {});
    if (error === COPY_WRITE_FAILED) return { copy: null, reason: 'write-failed' };
    throw error;
  }
};

/**
 * Where one entry's data waits while a Session Save writes its file (J3): a
 * new copy, so a large Session waits in the browser's storage, not in memory;
 * or no copy, because the browser lets none be written (`unavailable`: the
 * stream is unread, and the caller keeps it another way) or refused the write,
 * as when its storage is full (`write-failed`: what was written is deleted).
 * @typedef {{ copy: InputCopy } | { copy: null, reason: 'unavailable' | 'write-failed' }} SavedEntryData
 */
/**
 * @param {ReadableStream<Uint8Array>} stream
 * @returns {Promise<SavedEntryData>}
 */
export const keepSaveData = async (stream) => {
  const writer = await openCopyWriter();
  return writer ? writeStreamCopy(writer, stream) : { copy: null, reason: 'unavailable' };
};

/**
 * A stream that gives `source`'s bytes and calls `end` once it ends, fails or is cancelled.
 * @param {ReadableStream<Uint8Array>} source @param {() => void} end
 * @returns {ReadableStream<Uint8Array>}
 */
const streamUntil = (source, end) => {
  const reader = source.getReader();
  return new ReadableStream({
    async pull(controller) {
      const step = await reader.read().catch((error) => { end(); throw error; });
      if (!step.done) {
        controller.enqueue(step.value);
        return;
      }
      end();
      controller.close();
    },
    cancel(reason) {
      end();
      return reader.cancel(reason);
    }
  });
};

/**
 * A stream of a File's own bytes, when it has no copy. A File read before is
 * checked as it streams against the stored facts of its first read, so one
 * whose bytes changed on disk since rejects with `INPUT_UNREADABLE`
 * (PD-OI-121), as a refused read does; a File never read (no copy store) has
 * nothing to be checked against and streams as it is.
 * @param {object | null} key @param {ReadableFile} file
 * @returns {Promise<ReadableStream<Uint8Array>>}
 */
const checkedFileStream = async (key, file) => {
  const known = key ? fileContentFacts.get(key) : undefined;
  const expected = known ? await known : null;
  const changed = () => diagnosticError('INPUT_UNREADABLE', {}, { stage: 'read' });
  if (expected && typeof file.size === 'number' && file.size !== expected.byteLength) throw changed();
  const source = typeof file.stream === 'function' ? file.stream() : new Blob([/** @type {Uint8Array<ArrayBuffer>} */ (await readSource(file))]).stream();
  const reader = source.getReader();
  const hasher = expected ? pieceHasher() : null;
  return new ReadableStream({
    async pull(controller) {
      const step = await reader.read().catch((error) => { throw refusedReadError(error); });
      if (!step.done) {
        await hasher?.add(step.value);
        controller.enqueue(step.value);
        return;
      }
      const facts = await hasher?.finish();
      if (expected && facts && (facts.byteLength !== expected.byteLength || facts.token !== expected.token)) throw changed();
      controller.close();
    },
    cancel(reason) {
      return reader.cancel(reason);
    }
  });
};

/**
 * A stream of a File's bytes as chosen, which a Session file stores: its copy
 * (PD-OI-121), read only as it is streamed, or, with no copy or when the copy
 * is gone, the File itself, checked as it streams (`checkedFileStream`). The
 * copy is held until the stream ends.
 * @param {FileLike} file
 * @returns {Promise<ReadableStream<Uint8Array>>}
 */
const readStoredFileStream = async (file) => {
  const key = asCacheKey(file);
  const end = holdFileCopy(file);
  try {
    const copy = key ? await fileCopies.get(key) : null;
    /** @type {ReadableStream<Uint8Array> | null} */
    let source = null;
    if (copy) {
      try {
        source = (await copy.blob()).stream();
      } catch {
        // The copy is gone (the site's storage was cleared): the File is read itself.
      }
    }
    return streamUntil(source ?? await checkedFileStream(key, /** @type {ReadableFile} */ (file)), end);
  } catch (error) {
    end();
    throw error;
  }
};

/**
 * What a Session file stores for a File (`StoredResourceBytes`): its bytes as
 * chosen, a gzip file's as they are (receipt session.file-container), streamed
 * one piece at a time and never held whole, and staged as the File's content
 * Blob (`inputContentBlob`). Their size and token are known from the File's
 * first read; a File no read has read (no copy store) has its own size, and no
 * token until a stream of it is hashed.
 * @param {FileLike} file
 * @returns {Promise<StoredResourceBytes>}
 */
export const storedFileBytes = async (file) => {
  const key = asCacheKey(file);
  const readable = /** @type {ReadableFile} */ (file);
  const open = () => readStoredFileStream(file);
  const blob = () => inputContentBlob(file);
  const known = key ? fileContentFacts.get(key) : undefined;
  if (!known && typeof readable.stream === 'function' && typeof readable.size === 'number') {
    return { size: readable.size, open, token: null, blob };
  }
  const { token, byteLength } = await (known ?? readFileContentFacts(file));
  return { size: byteLength, open, token: async () => token, blob };
};

/**
 * The first read of a File: the one read of its own bytes, as chosen. For a
 * File the input store holds or a use holds, a copy of them (a gzip file's
 * gzip bytes) is written in the background (PD-OI-121). While it is written, a
 * gzip file's expanded size is counted against the cap, keeping nothing
 * expanded (`checkInputExpansion`); readers get the bytes once that count
 * ends, and their facts follow from them. A gzip file past the cap rejects
 * every reader, its facts and its copy with the refusal (PD-OI-116), and
 * keeps nothing, so its next read is a first read again.
 * @param {object} key @param {ReadableFile} file
 * @returns {{ bytes: Promise<Uint8Array>, done: Promise<unknown> }}
 */
const readFirst = (key, file) => {
  const chosen = readSource(file);
  const bytes = chosen.then(async (raw) => {
    await checkInputExpansion(String(file.name ?? ''), raw);
    return raw;
  });
  // The facts start as the check ends, before a reader of the bytes can transfer them.
  const facts = rememberFacts(key, bytes.then(factsOfBytes));
  const store = keptFiles.has(key) || copyUses.get(key) ? copyStore : null;
  // `writeCopy` closes the copy only under the facts, so a refused file's copy is deleted.
  const copy = store ? chosen.then((raw) => writeCopy(store, raw, facts), () => null) : Promise.resolve(null);
  fileCopies.set(key, copy);
  return { bytes, done: Promise.allSettled([bytes, facts, copy]) };
};

/**
 * Reads `file`'s bytes as chosen (cached under `key`): its first read, or a
 * later read through its copy, and shares them with the reads the current
 * holds cover (`holdReadBytes`).
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

// The reads of a Session resource's views are shared as a File's are: while a
// hold covers them, every view of the resource takes the same read, stored or
// plain (OV-2756), and the page drops it with the File reads once no hold is
// left. The backing keeps a read on its own only when its bytes rest in memory.
/** @type {WeakMap<object, { stored: object, plain: object }>} */
const viewReadKeys = new WeakMap();
/**
 * @param {SessionResourceFileView} view
 * @param {'stored' | 'plain'} kind
 * @param {(view: SessionResourceFileView) => Promise<Uint8Array> | null} read The backing's read; a view has a
 *   backing, so it is never null.
 * @returns {Promise<Uint8Array>}
 */
const sharedViewRead = (view, kind, read) => {
  const resource = sessionResourceSource(view)?.descriptor;
  if (!resource) return /** @type {Promise<Uint8Array>} */ (read(view));
  let keys = viewReadKeys.get(resource);
  if (!keys) viewReadKeys.set(resource, keys = { stored: {}, plain: {} });
  const key = keys[kind];
  const shared = sharedReads.get(key);
  if (shared) return shared;
  const pending = /** @type {Promise<Uint8Array>} */ (read(view));
  sharedReads.set(key, pending);
  const settle = () => {
    settledReads.add(pending);
    dropSharedReads();
  };
  pending.then(settle, settle);
  return pending;
};

/**
 * The stored bytes of a Session resource view (the size its descriptor
 * counts); a view has a backing, so its read is never null.
 * @param {SessionResourceFileView} view
 * @returns {Promise<Uint8Array>}
 */
const readStoredViewBytes = (view) => sharedViewRead(view, 'stored', readSessionResourceStoredBytes);

/**
 * The bytes of an uploaded File as chosen, or of a Session resource view as
 * stored, the one place they are read: what their descriptor counts, what
 * Generate and a Worker helper send (Python expands a gzip input), and what a
 * Session file stores. A read the browser refuses rejects with `INPUT_UNREADABLE`.
 * @param {unknown} file A `File`, a Session resource view, or any value (rejected unless it has `arrayBuffer()`).
 * @returns {Promise<Uint8Array>}
 */
export const readStoredBytes = (file) => {
  try {
    const key = asCacheKey(file);
    // A Session resource view caches its bytes in the backing it shares with
    // every other view of the same resource, so a transfer through one view
    // releases them for all (OV-304). The cache below is for ordinary Files.
    if (key && isSessionResourceFileView(file)) return readStoredViewBytes(file);
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
 * The plain bytes of an uploaded File or a Session resource view, which the
 * page's readers read: a gzip input (NCBI's `.gbff.gz`, `.fna.gz`, `.gff.gz`),
 * told by its first two bytes, expands at every read under the cap and is
 * never kept expanded (`expandInputBytes`); any other input is its bytes as
 * chosen. A read the browser refuses rejects with `INPUT_UNREADABLE`, as a
 * gzip file over the expanded cap does.
 * @param {unknown} file A `File`, a Session resource view, or any value (rejected unless it has `arrayBuffer()`).
 * @returns {Promise<Uint8Array>}
 */
export const readFileBytes = (file) => {
  if (isSessionResourceFileView(file)) return sharedViewRead(file, 'plain', readSessionResourceBytes);
  return readStoredBytes(file).then((bytes) => expandInputBytes(
    String(/** @type {{ readonly name?: unknown }} */ (file).name ?? ''), bytes));
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
 * The content facts of a File or Session resource view (`ContentFacts`),
 * computed once per object from its first read; known facts read nothing.
 * @param {FileLike} file
 * @returns {Promise<ContentFacts>}
 */
export const readFileContentFacts = (file) => {
  const key = asCacheKey(file);
  if (!key) return Promise.reject(new TypeError('A File-like object is required.'));
  const known = fileContentFacts.get(key);
  if (known) return known;
  // A Session resource view's facts are those of the stored bytes its descriptor counts: known when their
  // copy was hashed as it was written (a Load's, OV-2753), else hashed from a read.
  const stored = sessionResourceStoredBytes(file);
  if (stored) {
    const { size, token } = stored;
    return rememberFacts(key, token
      ? token().then((known) => ({ token: known, byteLength: size }))
      : readStoredViewBytes(/** @type {SessionResourceFileView} */ (file)).then(factsOfBytes));
  }
  // A File's first read keeps its facts.
  const bytes = readStoredBytes(file);
  return fileContentFacts.get(key) ?? rememberFacts(key, bytes.then(factsOfBytes));
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
  const waiting = holdWaiters.get(key);
  holdWaiters.delete(key);
  waiting?.forEach((resume) => resume(fileContentFacts.get(key)));
};

/**
 * The chosen Files whose refusal waits for the input store to hold them: a
 * field that asks before it binds the file (another drawing uses the file it
 * replaces; a stack row reaches the ring it replaces) holds it only after the
 * answer, and never after Cancel.
 * @type {WeakMap<object, ((facts: Promise<unknown> | undefined) => void)[]>}
 */
const holdWaiters = new WeakMap();

/**
 * The refusal of a chosen File's first read, once it ends: a gzip file whose
 * expanded size passes the cap (PD-OI-116), or null. It starts no read: the
 * input store reads a File when it first holds it (`keepFileCopy`), so a File
 * it does not hold yet is answered once it does, as after a question its field
 * asks before binding the file. With no copy store no read starts at the
 * choice, and it gives null.
 * @param {unknown} file
 * @returns {Promise<unknown>}
 */
export const chosenFileRefusal = (file) => {
  const key = asCacheKey(file);
  if (!key || isSessionResourceFileView(file)) return Promise.resolve(null);
  /** @param {Promise<unknown> | undefined} facts */
  const refusalOf = (facts) => (facts
    ? facts.then(() => null, (error) => (isInputExpansionRefusal(error) ? error : null))
    : Promise.resolve(null));
  const facts = fileContentFacts.get(key);
  if (facts || keptFiles.has(key) || !copyStore) return refusalOf(facts);
  if (typeof /** @type {{ readonly arrayBuffer?: unknown }} */ (key).arrayBuffer !== 'function') return Promise.resolve(null);
  return new Promise((resolve) => {
    holdWaiters.set(key, [...(holdWaiters.get(key) ?? []), (held) => resolve(refusalOf(held))]);
  });
};

/**
 * Shows a chosen file's read refusal, a gzip file past the expanded cap
 * (PD-OI-116), in the page's error panel; any other read failure is not its to
 * show (its field keeps the file, and Generate reports it). The one step every
 * field that refuses a chosen file takes.
 * @param {{ value: unknown }} errorLog The page's error panel (`state.errorLog`).
 * @param {unknown} error
 */
export const showReadRefusal = (errorLog, error) => {
  if (!isInputExpansionRefusal(error)) return;
  errorLog.value = normalizeUserFacingError(/** @type {Record<string, unknown>} */ (error), { stage: 'read' });
};

/**
 * Shows the refusal of a chosen File's first read (`chosenFileRefusal`) once
 * it ends (`showReadRefusal`).
 * @param {{ value: unknown }} errorLog The page's error panel (`state.errorLog`).
 * @param {unknown} file
 * @returns {Promise<void>}
 */
export const showChosenFileRefusal = (errorLog, file) => chosenFileRefusal(file).then((error) => showReadRefusal(errorLog, error));

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

/**
 * Keeps the copy of `file` while a state the input store does not see can
 * read it again (a failed import that offers Retry), even after the input
 * store lets it go. A File no other use holds gets a copy when this use holds
 * it during its first read.
 * @param {unknown} file
 * @returns {() => void} Ends this use; a second call does nothing.
 */
export const holdFileCopy = (file) => {
  const key = asCacheKey(file);
  if (!key) return () => {};
  copyUses.set(key, (copyUses.get(key) ?? 0) + 1);
  let held = true;
  return () => {
    if (!held) return;
    held = false;
    copyUses.set(key, (copyUses.get(key) ?? 1) - 1);
    if (fileCopies.has(key)) releasedFiles.add(key);
    dropSharedReads();
  };
};

/**
 * A copy a failed import keeps for its Retry: the error that offers Retry and
 * the drawing it failed in.
 * @typedef {{ readonly error: unknown, readonly drawing: object, readonly end: () => void }} RetryCopy
 */
/** @type {Set<RetryCopy>} */
const retryCopies = new Set();
/** @param {RetryCopy | null} held */
const endRetryCopy = (held) => {
  if (held && retryCopies.delete(held)) held.end();
};

/**
 * The copy that one kind of import's failure offers to Retry (PD-OI-121), so
 * Retry reads the content chosen even after the file changed: each call keeps
 * the copy of `file` for the failure `offer` names and lets go of the one kept
 * before; null (no failure offers Retry) lets it go. A copy also goes once its
 * Retry can no longer be offered (`endStaleRetryCopies`).
 * @returns {(file: unknown, offer?: { error: unknown, drawing: object }) => void}
 */
export const createRetryCopyHold = () => {
  /** @type {RetryCopy | null} */
  let held = null;
  return (file, offer) => {
    const next = file && offer ? { error: offer.error, drawing: offer.drawing, end: holdFileCopy(file) } : null;
    if (next) retryCopies.add(next);
    endRetryCopy(held);
    held = next;
  };
};

/**
 * Ends the Retry copies whose Retry can no longer be offered: their error is
 * no longer the one shown (cleared or replaced), or a Load or History restore
 * replaced their drawing.
 * @param {unknown} shownError @param {readonly object[]} drawings The drawings now.
 */
export const endStaleRetryCopies = (shownError, drawings) => {
  for (const held of [...retryCopies]) if (held.error !== shownError || !drawings.includes(held.drawing)) endRetryCopy(held);
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

/**
 * An input's bytes as chosen (a File's) or as stored (a Session resource
 * view's) as the diagram Worker reads them (OV-2749): a Blob, which
 * postMessage passes as a handle whatever its size, so no input's bytes cross
 * it. They are the bytes its descriptor counts and its content token hashes;
 * a gzip input's are its gzip bytes, which Python expands under its cap. A File
 * with a copy gives the copy file (PD-OI-121), and a view whose stored bytes
 * rest in a file gives their Blob (a Load's copy, OV-2753), neither of which
 * the page reads; any other input gives the bytes `readStoredBytes` reads (a
 * File without a copy checked against its first read) in a Blob.
 * @param {FileLike} file
 * @returns {Promise<Blob>}
 */
export const inputContentBlob = async (file) => {
  const stored = sessionResourceStoredBytes(file);
  if (stored?.blob) return stored.blob();
  const key = asCacheKey(file);
  // A first read writes the copy: whether the File has one is known once it ends.
  // Meanwhile a hold keeps the bytes that read shares, in case there is none.
  const releaseReads = holdReadBytes();
  try {
    const copy = key ? await fileCopies.get(key) : null;
    if (copy) {
      try {
        return await copy.blob();
      } catch {
        // The copy is gone (the site's storage was cleared): the File is read itself.
      }
    }
    return new Blob([/** @type {Uint8Array<ArrayBuffer>} */ (await readStoredBytes(file))]);
  } finally {
    releaseReads();
  }
};

/**
 * The bytes a History step keeps alive by naming `file`, as the History budget
 * counts them (OV-2750): none for a File, whose bytes stay on disk (the chosen
 * file, or its copy) and whose bytes this cache holds and lets go of on its
 * own, nor for a Session resource view whose stored bytes rest in a file, of
 * which the page keeps no read (OV-2754); the whole size of any other view, or
 * of any other object taken for a file, whose bytes are in memory.
 * @param {{ readonly size?: unknown } | null | undefined} file
 * @returns {number}
 */
export const residentBytesOf = (file) => {
  const key = asCacheKey(file);
  if (!key || key instanceof Blob || sessionResourceStoredBytes(key)?.blob) return 0;
  return Number(/** @type {{ readonly size?: unknown }} */ (key).size) || 0;
};

/**
 * A Worker helper's file for `role`: `file`'s content Blob (`inputContentBlob`)
 * and its name, which a helper's refusal of the file names (a gzip input over
 * its cap, PD-OI-116).
 * @param {string} role @param {FileLike} file
 * @returns {Promise<WorkerInputFile>}
 */
export const helperInputFile = async (role, file) => ({
  role, name: String(file.name ?? ''), blob: await inputContentBlob(file)
});
