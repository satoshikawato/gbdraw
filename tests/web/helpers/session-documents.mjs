// @ts-check
// Session documents as node tests write and read them. A test builds a
// Session as a JSON document whose inputs are base64 descriptors (as Sessions
// 27-50 held them); `sessionFileOf` turns a Session 51 or later into the
// Session file container a user loads (each input an entry, its descriptor the
// container's), and an older one into its JSON file. `readSavedSession` reads
// the file a Save wrote back into that JSON form, each entry's bytes as base64.
import { crc32, deflateRawSync } from 'node:zlib';
import { SESSION_DRAWINGS_LAYOUT } from '../../../gbdraw/web/js/contracts/schema.generated.js';
import { sessionContainerResources } from '../../../gbdraw/web/js/services/session-file.js';
import {
  isSessionContainerFile, openSessionContainerEntry, readSessionContainer, writeSessionContainer
} from '../../../gbdraw/web/js/services/session-container.js';

/** @typedef {Record<string, any>} SessionJson */

/**
 * The Session file of `document`: a container for a Session 51 or later,
 * else its JSON text.
 * @param {SessionJson} document
 * @param {(written: SessionJson) => string} [textOf] The text of the document as the file holds it (a test may write
 *   a number JSON.stringify cannot, such as `1e400`).
 * @returns {Promise<Blob>}
 */
export const sessionFileOf = async (document, textOf = JSON.stringify) => {
  if (!(typeof document?.version === 'number' && document.version >= SESSION_DRAWINGS_LAYOUT.version)) {
    return new Blob([textOf(document)], { type: 'application/json' });
  }
  /** @type {Record<string, Uint8Array>} */
  const bytes = {};
  /** @type {Record<string, any>} */
  const sized = {};
  // A document without a resources object, or with a resource entry that is no
  // descriptor (a case the reader refuses or repairs), cannot go through the
  // writer: its container is built by hand, those entries kept as they are.
  const declared = document.resources && typeof document.resources === 'object' && !Array.isArray(document.resources)
    ? document.resources : null;
  let damaged = declared === null;
  for (const [id, descriptor] of Object.entries(declared ?? {})) {
    if (!descriptor || typeof descriptor !== 'object' || Array.isArray(descriptor)) { damaged = true; continue; }
    bytes[id] = typeof descriptor.data === 'string' ? Buffer.from(descriptor.data, 'base64') : new Uint8Array();
    sized[id] = { ...descriptor, size: bytes[id].length };
  }
  if (damaged) {
    const written = declared ? { ...document, resources: { ...declared, ...sessionContainerResources(sized) } } : document;
    return handBuiltContainer(textOf(written), bytes);
  }
  const written = { ...document, resources: sessionContainerResources(sized) };
  return writeSessionContainer({
    document: written,
    documentBytes: () => new Blob([textOf(written)]).stream(),
    readResource: async (id) => new Blob([bytes[id]]).stream()
  });
};

/**
 * A container of the document text `documentText` with an entry for each of
 * `entries`, written without the checks of `writeSessionContainer`: the
 * document entry first, every entry deflated.
 * @param {string} documentText
 * @param {Record<string, Uint8Array>} entries
 * @returns {Blob}
 */
const handBuiltContainer = (documentText, entries) => {
  const files = [['gbdraw-session.json', Buffer.from(documentText)],
    ...Object.entries(entries).map(([id, data]) => [`resources/${id}`, Buffer.from(data)])];
  /** @type {Buffer[]} */
  const parts = [];
  /** @type {Buffer[]} */
  const central = [];
  let offset = 0;
  for (const [name, data] of files) {
    const nameBytes = Buffer.from(String(name));
    const packed = deflateRawSync(data);
    const crc = crc32(data);
    const local = Buffer.alloc(30);
    local.writeUInt32LE(0x04034b50, 0); local.writeUInt16LE(20, 4); local.writeUInt16LE(0, 6); local.writeUInt16LE(8, 8);
    local.writeUInt16LE(0, 10); local.writeUInt16LE(0x21, 12); local.writeUInt32LE(crc, 14);
    local.writeUInt32LE(packed.length, 18); local.writeUInt32LE(data.length, 22); local.writeUInt16LE(nameBytes.length, 26);
    const entry = Buffer.alloc(46);
    entry.writeUInt32LE(0x02014b50, 0); entry.writeUInt16LE(20, 4); entry.writeUInt16LE(20, 6); entry.writeUInt16LE(0, 8);
    entry.writeUInt16LE(8, 10); entry.writeUInt16LE(0, 12); entry.writeUInt16LE(0x21, 14); entry.writeUInt32LE(crc, 16);
    entry.writeUInt32LE(packed.length, 20); entry.writeUInt32LE(data.length, 24); entry.writeUInt16LE(nameBytes.length, 28);
    entry.writeUInt32LE((0o100644 << 16) >>> 0, 38); entry.writeUInt32LE(offset, 42);
    parts.push(local, nameBytes, packed);
    central.push(entry, nameBytes);
    offset += 30 + nameBytes.length + packed.length;
  }
  const directory = Buffer.concat(central);
  const end = Buffer.alloc(22);
  end.writeUInt32LE(0x06054b50, 0); end.writeUInt16LE(files.length, 8); end.writeUInt16LE(files.length, 10);
  end.writeUInt32LE(directory.length, 12); end.writeUInt32LE(offset, 16);
  return new Blob([...parts, directory, end]);
};

/**
 * The Load event of `document`'s Session file, as a file picker gives it.
 * @param {SessionJson} document
 */
export const sessionLoadEvent = async (document) => ({
  target: { files: [await sessionFileOf(document)], value: 'selected' }
});

/**
 * The document a Session file holds, each container entry back as a base64
 * descriptor; a JSON file (plain or gzip) is parsed as it is.
 * @param {Blob} blob
 * @returns {Promise<SessionJson>}
 */
export const readSavedSession = async (blob) => {
  if (!await isSessionContainerFile(blob)) {
    const { gunzipSync } = await import('node:zlib');
    const raw = Buffer.from(await blob.arrayBuffer());
    return JSON.parse((raw[0] === 0x1f && raw[1] === 0x8b ? gunzipSync(raw) : raw).toString('utf8'));
  }
  const read = await readSessionContainer(blob);
  const document = /** @type {SessionJson} */ (read.document);
  for (const [id, entry] of Object.entries(read.resources)) {
    const data = Buffer.from(await new Response(openSessionContainerEntry(blob, entry)).arrayBuffer());
    document.resources[id] = { ...document.resources[id], encoding: 'base64', data: data.toString('base64') };
  }
  return document;
};

/**
 * The Gallery Session `id` (`<id>.gbdraw`) in the JSON form.
 * @param {string} id
 * @returns {Promise<SessionJson>}
 */
export const gallerySessionDocument = async (id) => {
  const { readFileSync } = await import('node:fs');
  const url = new URL(`../../../gbdraw/web/gallery/sessions/${id}.gbdraw`, import.meta.url);
  return readSavedSession(new Blob([readFileSync(url)]));
};

/**
 * The Load event of Session bytes a test holds (a container, or a JSON file
 * plain or gzip): a JSON Session 51 or later becomes its container; anything
 * else is handed over as it is.
 * @param {Uint8Array | string} bytes
 */
export const sessionLoadEventOfBytes = async (bytes) => {
  const raw = typeof bytes === 'string' ? Buffer.from(bytes) : Buffer.from(bytes);
  if (raw[0] === 0x50 && raw[1] === 0x4b) return { target: { files: [new Blob([raw])], value: 'selected' } };
  try {
    const document = await readSavedSession(new Blob([raw]));
    return sessionLoadEvent(document);
  } catch {
    return { target: { files: [new Blob([raw], { type: 'application/json' })], value: 'selected' } };
  }
};

/**
 * The JSON text of the Session a file holds (container, JSON or gzip JSON).
 * @param {Uint8Array} bytes
 * @returns {Promise<Buffer>}
 */
export const sessionJsonBytes = async (bytes) => Buffer.from(JSON.stringify(await readSavedSession(new Blob([bytes]))));
