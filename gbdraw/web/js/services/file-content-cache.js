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
// them: until the copy is written, or for good when it has none. They stay
// until a transfer takes them: expanding again would be slower than keeping them.
/** @type {WeakSet<object>} */
const expandedFiles = new WeakSet();
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
 * @returns {Promise<FileContentFacts>}
 */
const factsOfBytes = async (bytes) => {
  const byteLength = bytes.byteLength;
  return { token: await contentTokenOfBytes(bytes), byteLength };
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
 * The plain bytes of a File's bytes, the one setup of an input's expansion: a
 * gzip file expands to at most Python's `MAX_EXPANDED_INPUT_BYTES`, and past it
 * rejects with `INPUT_UNREADABLE` (`GZIP_EXPANDED_SIZE`), naming the file and
 * the cap in MiB.
 * @param {string} name The File's name.
 * @param {Uint8Array<ArrayBuffer>} bytes
 * @returns {Promise<Uint8Array>}
 */
const expandInputBytes = (name, bytes) => {
  const maxBytes = testByteCap('maxExpandedInputBytes') ?? MAX_EXPANDED_INPUT_BYTES;
  return expandGzipBytes(bytes, maxBytes, () => diagnosticError('INPUT_UNREADABLE',
    { reason: 'GZIP_EXPANDED_SIZE', inputName: name, limit: Math.floor(maxBytes / MIB) }, { stage: 'read' }));
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
 * Writes `plain` to a new copy in 16 MiB pieces, one at a time, while its facts
 * are hashed, and closes it under their token; null when no copy could be
 * written, or when it holds fewer bytes than the facts count.
 * Every piece goes through one buffer, which the Worker gives back, so the copy
 * holds one piece beside `plain`.
 * @param {InputCopyStore} store @param {Uint8Array} plain @param {Promise<FileContentFacts>} facts
 * @returns {Promise<InputCopy | null>}
 */
const writeCopy = async (store, plain, facts) => {
  const writer = await store.openWriter().catch(() => null);
  if (!writer) return null;
  try {
    /** @type {ArrayBuffer | null} */
    let spare = null;
    let written = 0;
    for (let start = 0; start === 0 || start < plain.byteLength; start += TOKEN_PIECE_BYTES) {
      const piece = plain.subarray(start, start + TOKEN_PIECE_BYTES);
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
 * @typedef {{ copy: InputCopy, facts: FileContentFacts }
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
  /** @type {ArrayBuffer[]} */
  const digests = [];
  let byteLength = 0;
  let buffer = new ArrayBuffer(TOKEN_PIECE_BYTES);
  let filled = 0;
  // Web Crypto copies a piece when the digest starts, so the write may transfer it at once.
  const writePiece = async () => {
    const piece = filled === buffer.byteLength ? buffer : buffer.slice(0, filled);
    const [digest, given] = await Promise.all([
      sha256Digest(new Uint8Array(piece)),
      writer.write(piece).catch(() => { throw COPY_WRITE_FAILED; })
    ]);
    digests.push(digest);
    byteLength += filled;
    if (piece === buffer) buffer = given.byteLength === TOKEN_PIECE_BYTES ? given : new ArrayBuffer(TOKEN_PIECE_BYTES);
    filled = 0;
  };
  try {
    for (let step = await reader.read(); !step.done; step = await reader.read()) {
      for (let chunk = step.value; chunk.byteLength > 0;) {
        const taken = Math.min(chunk.byteLength, buffer.byteLength - filled);
        new Uint8Array(buffer, filled, taken).set(chunk.subarray(0, taken));
        filled += taken;
        chunk = chunk.subarray(taken);
        if (filled === buffer.byteLength) await writePiece();
      }
    }
    // The last piece, or the one empty piece of empty content (`contentTokenOfBytes`).
    if (filled > 0 || digests.length === 0) await writePiece();
    const token = await tokenOfPieceDigests(digests);
    const copy = await writer.finish(token).catch(() => { throw COPY_WRITE_FAILED; });
    return { copy, facts: { token, byteLength } };
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
 * A stream of the plain bytes of a File or Session resource view (`readFileBlob`)
 * that a Session Save writes; the File's copy is held until the stream ends.
 * @param {FileLike} file
 * @returns {Promise<ReadableStream<Uint8Array>>}
 */
export const readFileStream = async (file) => {
  const end = holdFileCopy(file);
  try {
    const reader = (await readFileBlob(file)).stream().getReader();
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
  } catch (error) {
    end();
    throw error;
  }
};

/**
 * The plain bytes of a File or Session resource view as a Blob that a Session
 * Save streams: the File's copy (PD-OI-121), read only as it is streamed, or,
 * with no copy or when the copy is gone, the bytes `readFileBytes` reads. A
 * caller that streams it while the File may be let go holds the copy
 * (`holdFileCopy`).
 * @param {FileLike} file
 * @returns {Promise<Blob>}
 */
export const readFileBlob = async (file) => {
  const key = asCacheKey(file);
  const copy = key ? await fileCopies.get(key) : null;
  if (copy) {
    try {
      return await copy.blob();
    } catch {
      // The copy is gone (the site's storage was cleared): the File is read itself.
    }
  }
  return new Blob([/** @type {Uint8Array<ArrayBuffer>} */ (await readFileBytes(file))]);
};

/**
 * The first read of a File: the one read of its own bytes. Readers get them, a
 * gzip file's expanded, as soon as they are read, and their facts follow from
 * them. For a File the input store holds or a use holds, a copy is then written from them in
 * the background (PD-OI-121); a gzip file's expanded bytes are dropped once it
 * is written, and kept when no copy could be.
 * @param {object} key @param {ReadableFile} file
 * @returns {{ bytes: Promise<Uint8Array>, done: Promise<unknown> }}
 */
const readFirst = (key, file) => {
  const bytes = readSource(file).then(async (raw) => {
    const plain = await expandInputBytes(String(file.name ?? ''), raw);
    if (plain !== raw) expandedFiles.add(key);
    return plain;
  });
  const facts = rememberFacts(key, bytes.then(factsOfBytes));
  const store = keptFiles.has(key) || copyUses.get(key) ? copyStore : null;
  const copy = store ? bytes.then((plain) => writeCopy(store, plain, facts), () => null) : Promise.resolve(null);
  fileCopies.set(key, copy);
  const copied = copy.then((written) => {
    if (written) expandedFiles.delete(key);
  });
  return { bytes, done: Promise.allSettled([bytes, facts, copied]) };
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
  // A first read hashes the bytes it read and copies them; a transfer detaches them once both end.
  const key = asCacheKey(file);
  if (key) await Promise.allSettled([fileContentFacts.get(key), fileCopies.get(key)]);
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
