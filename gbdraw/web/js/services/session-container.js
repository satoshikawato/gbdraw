// @ts-check
// The Session file container (`<title>.gbdraw`): one ZIP whose entry
// `gbdraw-session.json` holds the settings and results, and whose entries
// `resources/<id>` hold one input file each. Python owns the format
// (`gbdraw/session_container.py`); this reader and writer are its registered
// twin (`tools/twins.json`, R9), pinned to it by the shared vectors
// (`tests/fixtures/session_container_vectors.json`). Entries stream: no input
// is held whole, and none becomes a base64 text.
//
// The reader reads the central directory from the end of the file, then the
// document, and checks every entry's size and CRC-32 while it streams. It
// refuses, as damaged, a file whose directory is missing, cut off or
// unreadable, two entries of one name, an encrypted entry or one compressed
// other than stored or deflated, a document or resource whose size or CRC-32
// does not match, and a resource the document declares with no entry or with
// another size; a ZIP without `gbdraw-session.json` is refused as no Session,
// and a document that is not UTF-8 JSON as not valid JSON. An entry name is
// only compared with `gbdraw-session.json` and `resources/<id>`: it never
// becomes a path, a resource is found by its id, and entries the document does
// not declare are not read.
import { crc32Update, zipCentralRecord, zipEndRecords, zipLocalHeader } from '../utils/zip.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { SESSION_DOCUMENT_ENTRY, SESSION_RESOURCE_ENTRY_PREFIX } from '../contracts/schema.generated.js';
import { GZIP_SIGNATURE } from './byte-utils.js';
import { assertUniqueJsonKeys } from './safe-object-keys.js';

/** @import { ZipEntryFields } from '../utils/zip.js' */

const STORED = 0;
const DEFLATED = 8;
const ZIP32_MARK = 0xFFFFFFFF;
const END_RECORD_BYTES = 22;
const MAX_COMMENT_BYTES = 0xFFFF;
const ZIP64_END_BYTES = 56;
const ZIP64_LOCATOR_BYTES = 20;
// The DOS time of every entry: 1980-01-01 00:00, so equal content writes equal bytes.
const DOS_TIME = 0;
const DOS_DATE = (1 << 5) | 1;

/**
 * Where one entry of a container is and what its directory record declares:
 * plain data, so a Worker that read the directory can hand it to the page,
 * which opens the entry from the same file (`openSessionContainerEntry`).
 * @typedef {Readonly<{
 *   name: string, method: number, crc: number, compressedSize: number, size: number, offset: number, end: number
 * }>} SessionContainerEntry
 */
/**
 * A container read by `readSessionContainer`: its parsed document (not yet
 * checked as a Session) and the entries of the document and of each resource
 * it declares, by resource id.
 * @typedef {object} SessionContainer
 * @property {unknown} document The JSON value of `gbdraw-session.json`.
 * @property {SessionContainerEntry} documentEntry
 * @property {Readonly<Record<string, SessionContainerEntry>>} resources The entry of each declared resource,
 *   whose size is the size its descriptor declares.
 * @property {number} inputBytes The declared sizes of every entry but `gbdraw-session.json`, added up: the
 *   browser storage a Load of its inputs needs (the receipt's "inputs need {need} MB").
 */
/**
 * Called with the size the directory declares for `gbdraw-session.json`
 * before it is read; it throws to refuse the file (the Web's 512 MiB limit).
 * @typedef {(documentBytes: number) => void} SessionDocumentAdmission
 */
/**
 * Keeps the bytes of one entry's data while the container is written, and
 * gives them back as a Blob (by default in memory; a caller may keep them in
 * the browser's storage instead).
 * @typedef {(data: ReadableStream<Uint8Array>) => Promise<Blob>} SessionContainerEntryKeeper
 */
/**
 * What a container is written from: the Session document (its `resources`
 * name the entries and their order, each with its `size`), the document's
 * JSON as UTF-8 bytes, and each resource's bytes by id.
 * @typedef {{
 *   document: { readonly resources: Readonly<Record<string, unknown>> },
 *   documentBytes: () => ReadableStream<Uint8Array>,
 *   readResource: (resourceId: string) => ReadableStream<Uint8Array> | Promise<ReadableStream<Uint8Array>>
 * }} SessionContainerContent
 */
/** @typedef {{ -readonly [K in keyof SessionContainerEntry]: SessionContainerEntry[K] }} DirectoryEntry */

const invalidResource = () => diagnosticError('RESOURCE_INVALID', { reason: 'SIZE' }, { stage: 'write', operation: 'session-save' });

const damaged = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'CONTAINER_DAMAGED' },
  { stage: 'read', operation: 'session-load' });

/**
 * @param {Blob} blob @param {number} start @param {number} end
 * @returns {Promise<Uint8Array<ArrayBuffer>>}
 */
const bytesOf = async (blob, start, end) => new Uint8Array(await blob.slice(start, end).arrayBuffer());

/** @param {Uint8Array<ArrayBuffer>} bytes */
const viewOf = (bytes) => new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);

/** @param {DataView} view @param {number} at */
const bigAt = (view, at) => {
  const value = view.getBigUint64(at, true);
  if (value > BigInt(Number.MAX_SAFE_INTEGER)) throw damaged();
  return Number(value);
};

const utf8Names = new TextDecoder('utf-8', { fatal: true });
// A name without the UTF-8 flag: one character per byte, so two names are
// equal only when their bytes are; such a name never equals an ASCII entry
// name unless its bytes are that name.
const byteNames = new TextDecoder('latin1');

/**
 * Where the central directory is, from the end records; null when the file
 * has no end record.
 * @param {Blob} blob
 * @returns {Promise<{ count: number, size: number, offset: number, shift: number }>}
 */
const directoryOf = async (blob) => {
  const tailStart = Math.max(0, blob.size - END_RECORD_BYTES - MAX_COMMENT_BYTES);
  const tail = await bytesOf(blob, tailStart, blob.size);
  const view = viewOf(tail);
  let at = tail.length - END_RECORD_BYTES;
  for (; at >= 0; at -= 1) {
    if (view.getUint32(at, true) === 0x06054B50 && at + END_RECORD_BYTES + view.getUint16(at + 20, true) <= tail.length) break;
  }
  if (at < 0) throw damaged();
  if (view.getUint16(at + 4, true) !== 0 || view.getUint16(at + 6, true) !== 0) throw damaged();
  const endAt = tailStart + at;
  const locatorAt = endAt - ZIP64_LOCATOR_BYTES;
  const locator = locatorAt >= 0 ? await bytesOf(blob, locatorAt, endAt) : new Uint8Array(0);
  if (locator.length === ZIP64_LOCATOR_BYTES && viewOf(locator).getUint32(0, true) === 0x07064B50) {
    const locatorView = viewOf(locator);
    if (locatorView.getUint32(4, true) !== 0 || locatorView.getUint32(16, true) > 1) throw damaged();
    // The ZIP64 end record sits just before its locator, as Python's `zipfile` reads it.
    const recordAt = locatorAt - ZIP64_END_BYTES;
    if (recordAt < 0) throw damaged();
    const record = viewOf(await bytesOf(blob, recordAt, locatorAt));
    if (record.getUint32(0, true) !== 0x06064B50 || record.getUint32(16, true) !== 0 || record.getUint32(20, true) !== 0) {
      throw damaged();
    }
    const count = bigAt(record, 32);
    const size = bigAt(record, 40);
    const offset = bigAt(record, 48);
    return { count, size, offset, shift: recordAt - size - offset };
  }
  const count = view.getUint16(at + 10, true);
  const size = view.getUint32(at + 12, true);
  const offset = view.getUint32(at + 16, true);
  if (size === ZIP32_MARK || offset === ZIP32_MARK) throw damaged();
  return { count, size, offset, shift: endAt - size - offset };
};

/**
 * The entries of the central directory by name, each with the offset where the
 * next entry or the directory starts (an entry's data must end before it).
 * @param {Blob} blob
 * @returns {Promise<Map<string, DirectoryEntry>>}
 */
const entriesOf = async (blob) => {
  const { count, size, offset, shift } = await directoryOf(blob);
  // A file with bytes before the ZIP (`shift`) is read as Python's `zipfile` reads it.
  const start = offset + shift;
  if (shift < 0 || start < 0) throw damaged();
  const directory = await bytesOf(blob, start, start + size);
  if (directory.length !== size) throw damaged();
  const view = viewOf(directory);
  /** @type {Map<string, DirectoryEntry>} */
  const entries = new Map();
  let at = 0;
  while (at < size) {
    if (at + 46 > size || view.getUint32(at, true) !== 0x02014B50) throw damaged();
    const flags = view.getUint16(at + 8, true);
    const method = view.getUint16(at + 10, true);
    const nameLength = view.getUint16(at + 28, true);
    const extraLength = view.getUint16(at + 30, true);
    const commentLength = view.getUint16(at + 32, true);
    const next = at + 46 + nameLength + extraLength + commentLength;
    if (next > size) throw damaged();
    // Bit 0: encrypted; bit 6: strong encryption.
    if (flags & 0x41 || (method !== STORED && method !== DEFLATED)) throw damaged();
    const nameBytes = directory.slice(at + 46, at + 46 + nameLength);
    let name;
    try {
      name = (flags & 0x0800 ? utf8Names : byteNames).decode(nameBytes);
    } catch {
      throw damaged();
    }
    let compressedSize = view.getUint32(at + 20, true);
    let entrySize = view.getUint32(at + 24, true);
    let entryOffset = view.getUint32(at + 42, true);
    // The ZIP64 extra field holds, in this order, each of the three that is 0xFFFFFFFF above.
    const wanted = [entrySize, compressedSize, entryOffset].filter((value) => value === ZIP32_MARK).length;
    if (wanted) {
      let field = at + 46 + nameLength;
      const extraEnd = field + extraLength;
      for (; field + 4 <= extraEnd && view.getUint16(field, true) !== 1; field += 4 + view.getUint16(field + 2, true));
      if (field + 4 > extraEnd || view.getUint16(field + 2, true) < 8 * wanted || field + 4 + 8 * wanted > extraEnd) throw damaged();
      let value = field + 4;
      if (entrySize === ZIP32_MARK) { entrySize = bigAt(view, value); value += 8; }
      if (compressedSize === ZIP32_MARK) { compressedSize = bigAt(view, value); value += 8; }
      if (entryOffset === ZIP32_MARK) entryOffset = bigAt(view, value);
    }
    if (entries.has(name)) throw damaged();
    entries.set(name, { name, method, crc: view.getUint32(at + 16, true), compressedSize, size: entrySize, offset: entryOffset + shift, end: start });
    at = next;
  }
  if (entries.size !== count) throw damaged();
  const byOffset = [...entries.values()].sort((left, right) => left.offset - right.offset);
  byOffset.forEach((entry, index) => { entry.end = byOffset[index + 1]?.offset ?? start; });
  return entries;
};

/** @param {Uint8Array} left @param {Uint8Array} right */
const sameBytes = (left, right) => left.length === right.length && left.every((byte, index) => byte === right[index]);

/**
 * The content of `entry`, inflated when deflated and checked against the
 * directory's size and CRC-32 as it streams.
 * @param {Blob} blob @param {SessionContainerEntry} entry
 * @returns {ReadableStream<Uint8Array>}
 */
export const openSessionContainerEntry = (blob, entry) => {
  const nameBytes = new TextEncoder().encode(entry.name);
  /** @type {ReadableStreamDefaultReader<Uint8Array> | null} */
  let reader = null;
  let crc = 0;
  let size = 0;
  return new ReadableStream({
    async start() {
      const header = await bytesOf(blob, entry.offset, entry.offset + 30);
      if (header.length !== 30 || viewOf(header).getUint32(0, true) !== 0x04034B50) throw damaged();
      const nameLength = viewOf(header).getUint16(26, true);
      const dataStart = entry.offset + 30 + nameLength + viewOf(header).getUint16(28, true);
      if (!sameBytes(await bytesOf(blob, entry.offset + 30, entry.offset + 30 + nameLength), nameBytes)
        || dataStart + entry.compressedSize > entry.end) throw damaged();
      const data = blob.slice(dataStart, dataStart + entry.compressedSize).stream();
      // Python's `zipfile` reads a deflated entry of no bytes as empty content.
      reader = (entry.method === DEFLATED && entry.compressedSize > 0
        ? data.pipeThrough(new DecompressionStream('deflate-raw')) : data).getReader();
    },
    async pull(controller) {
      const source = /** @type {ReadableStreamDefaultReader<Uint8Array>} */ (reader);
      /** @type {ReadableStreamReadResult<Uint8Array>} */
      let step;
      try {
        step = await source.read();
      } catch {
        throw damaged();
      }
      if (step.done) {
        if (size !== entry.size || crc !== entry.crc) throw damaged();
        controller.close();
        return;
      }
      size += step.value.byteLength;
      if (size > entry.size) {
        source.cancel().catch(() => {});
        throw damaged();
      }
      crc = crc32Update(crc, step.value);
      controller.enqueue(step.value);
    },
    cancel(reason) {
      return reader?.cancel(reason);
    }
  });
};

/**
 * Whether `blob` is a container: content decides, as Python's reader decides
 * (a local header or an empty archive's end record first); any other file is a
 * JSON Session, plain or gzip.
 * @param {Blob} blob
 * @returns {Promise<boolean>}
 */
export const isSessionContainerFile = async (blob) => {
  if (blob.size < 4) return false;
  const signature = viewOf(await bytesOf(blob, 0, 4)).getUint32(0, true);
  return signature === 0x04034B50 || signature === 0x06054B50;
};

const refusedJson = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'JSON_FORMAT' },
  { stage: 'parse', operation: 'session-load' });

/**
 * Reads a container as Python's `read_session_container` does: its directory,
 * then `gbdraw-session.json` in full, then the entry of each resource the
 * document declares, checked against its declared size. A resource's own
 * bytes are checked when it is opened (`openSessionContainerEntry`).
 * @param {Blob} blob
 * @param {{ admitDocument?: SessionDocumentAdmission }} [options]
 * @returns {Promise<SessionContainer>}
 */
export const readSessionContainer = async (blob, { admitDocument = () => {} } = {}) => {
  const entries = await entriesOf(blob);
  const documentEntry = entries.get(SESSION_DOCUMENT_ENTRY);
  if (!documentEntry) {
    throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'CONTAINER_NO_DOCUMENT' },
      { stage: 'read', operation: 'session-load' });
  }
  admitDocument(documentEntry.size);
  const documentBytes = await new Response(openSessionContainerEntry(blob, documentEntry)).arrayBuffer();
  /** @type {unknown} */
  let document;
  try {
    const text = new TextDecoder('utf-8', { fatal: true }).decode(documentBytes);
    document = JSON.parse(text);
    // Python refuses a document with one key twice as invalid JSON (OV-2694).
    assertUniqueJsonKeys(text);
  } catch {
    throw refusedJson();
  }
  const declared = document && typeof document === 'object' && !Array.isArray(document)
    ? /** @type {{ resources?: unknown }} */ (document).resources : undefined;
  /** @type {Record<string, SessionContainerEntry>} */
  const resources = {};
  // A `resources` that is not an object declares none here; the Session reader refuses it.
  if (declared && typeof declared === 'object' && !Array.isArray(declared)) {
    for (const [resourceId, descriptor] of Object.entries(declared)) {
      const entry = entries.get(`${SESSION_RESOURCE_ENTRY_PREFIX}${resourceId}`);
      const size = descriptor && typeof descriptor === 'object' ? /** @type {{ size?: unknown }} */ (descriptor).size : undefined;
      if (!entry || entry === documentEntry || entry.size !== size) throw damaged();
      resources[resourceId] = entry;
    }
  }
  let inputBytes = 0;
  for (const entry of entries.values()) if (entry !== documentEntry) inputBytes += entry.size;
  return { document, documentEntry, resources, inputBytes };
};

/** @type {SessionContainerEntryKeeper} */
const keepInMemory = (data) => new Response(data).blob();

/**
 * The stream of `source`'s bytes, and whether they start with the gzip
 * signature (a gzip input is stored as it is): the first chunks are read to
 * tell, and given again at the start of the stream.
 * @param {ReadableStream<Uint8Array>} source
 * @returns {Promise<{ gzip: boolean, stream: ReadableStream<Uint8Array> }>}
 */
const withGzipSignature = async (source) => {
  const reader = source.getReader();
  /** @type {Uint8Array[]} */
  const head = [];
  let headBytes = 0;
  while (headBytes < GZIP_SIGNATURE.length) {
    const step = await reader.read();
    if (step.done) break;
    head.push(step.value);
    headBytes += step.value.byteLength;
  }
  const first = new Uint8Array(Math.min(headBytes, GZIP_SIGNATURE.length));
  for (let at = 0, index = 0; at < first.length; index += 1) {
    const part = head[index].subarray(0, first.length - at);
    first.set(part, at);
    at += part.length;
  }
  const gzip = first.length === GZIP_SIGNATURE.length && GZIP_SIGNATURE.every((byte, index) => first[index] === byte);
  const stream = new ReadableStream({
    async pull(controller) {
      const next = head.shift();
      if (next) {
        controller.enqueue(next);
        return;
      }
      const step = await reader.read();
      if (step.done) controller.close();
      else controller.enqueue(step.value);
    },
    cancel(reason) {
      return reader.cancel(reason);
    }
  });
  return { gzip, stream };
};

/**
 * Writes a container as Python's `write_session_container` does:
 * `gbdraw-session.json` first, then `resources/<id>` for each key of the
 * document's `resources`, in that order. Bytes that start with the gzip
 * signature are stored as they are; others, the document and an empty
 * resource are deflated. Each entry is read once: its size and CRC-32 are
 * taken as it streams to `keep`, and its header is written after. A resource
 * whose bytes are not the `size` its descriptor declares is refused.
 * @param {SessionContainerContent} content
 * @param {{ keep?: SessionContainerEntryKeeper, zip64?: boolean }} [options] `zip64` writes every size in the
 *   ZIP64 fields even when small (a test of the readers).
 * @returns {Promise<Blob>}
 */
export const writeSessionContainer = async ({ document, documentBytes, readResource }, { keep = keepInMemory, zip64 = false } = {}) => {
  /** @type {{ name: string, size: number | null, read: () => ReadableStream<Uint8Array> | Promise<ReadableStream<Uint8Array>> }[]} */
  const sources = [{ name: SESSION_DOCUMENT_ENTRY, size: null, read: documentBytes }];
  for (const [resourceId, descriptor] of Object.entries(document.resources)) {
    const size = descriptor && typeof descriptor === 'object' ? /** @type {{ size?: unknown }} */ (descriptor).size : undefined;
    if (typeof size !== 'number' || !Number.isSafeInteger(size) || size < 0) throw invalidResource();
    sources.push({ name: `${SESSION_RESOURCE_ENTRY_PREFIX}${resourceId}`, size, read: () => readResource(resourceId) });
  }
  const encoder = new TextEncoder();
  /** @type {BlobPart[]} */
  const parts = [];
  /** @type {ZipEntryFields[]} */
  const written = [];
  let offset = 0;
  for (const source of sources) {
    const { gzip, stream } = await withGzipSignature(await source.read());
    const stored = gzip && source.size !== null;
    let crc = 0;
    let size = 0;
    const measured = stream.pipeThrough(new TransformStream({
      transform(chunk, controller) {
        crc = crc32Update(crc, chunk);
        size += chunk.byteLength;
        controller.enqueue(chunk);
      }
    }));
    const data = await keep(stored ? measured : measured.pipeThrough(new CompressionStream('deflate-raw')));
    if (source.size !== null && size !== source.size) throw invalidResource();
    /** @type {ZipEntryFields} */
    const entry = {
      name: encoder.encode(source.name), method: stored ? STORED : DEFLATED, crc, compressedSize: data.size, size,
      offset, time: DOS_TIME, date: DOS_DATE, zip64
    };
    const header = zipLocalHeader(entry);
    parts.push(header, data);
    written.push(entry);
    offset += header.length + data.size;
  }
  const directory = written.map(zipCentralRecord);
  const directorySize = directory.reduce((total, record) => total + record.length, 0);
  parts.push(...directory, zipEndRecords({ count: written.length, directoryOffset: offset, directorySize, zip64 }));
  return new Blob(parts, { type: 'application/zip' });
};
