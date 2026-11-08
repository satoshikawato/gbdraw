// @ts-check
// Session documents as node tests write and read them. A test builds a
// Session as a JSON document whose inputs are base64 descriptors (as Sessions
// 27-50 held them); `sessionFileOf` turns a Session 51 or later into the
// Session file container a user loads (each input an entry, its descriptor the
// container's), and an older one into its JSON file. `readSavedSession` reads
// the file a Save wrote back into that JSON form, each entry's bytes as base64.
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
 * @returns {Promise<Blob>}
 */
export const sessionFileOf = async (document) => {
  if (!(typeof document?.version === 'number' && document.version >= SESSION_DRAWINGS_LAYOUT.version)) {
    return new Blob([JSON.stringify(document)], { type: 'application/json' });
  }
  /** @type {Record<string, Uint8Array>} */
  const bytes = {};
  /** @type {Record<string, any>} */
  const sized = {};
  // A document without a resources object (a case the reader refuses) is written as it is.
  const declared = document.resources && typeof document.resources === 'object' && !Array.isArray(document.resources)
    ? document.resources : null;
  for (const [id, descriptor] of Object.entries(declared ?? {})) {
    bytes[id] = typeof descriptor?.data === 'string' ? Buffer.from(descriptor.data, 'base64') : new Uint8Array();
    sized[id] = { ...descriptor, size: bytes[id].length };
  }
  const written = declared ? { ...document, resources: sessionContainerResources(sized) } : document;
  return writeSessionContainer({
    document: written,
    documentBytes: () => new Blob([JSON.stringify(written)]).stream(),
    readResource: async (id) => new Blob([bytes[id]]).stream()
  });
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
 * The Gallery Session `id` in the JSON form: its `<id>.gbdraw` when the
 * Gallery holds one, else its older JSON file.
 * @param {string} id
 * @returns {Promise<SessionJson>}
 */
export const gallerySessionDocument = async (id) => {
  const { existsSync, readFileSync } = await import('node:fs');
  const sessions = new URL('../../../gbdraw/web/gallery/sessions/', import.meta.url);
  for (const name of [`${id}.gbdraw`, `${id}.gbdraw-session.json`, `${id}.gbdraw-session.json.gz`]) {
    const url = new URL(name, sessions);
    if (existsSync(url)) return readSavedSession(new Blob([readFileSync(url)]));
  }
  throw new Error(`the Gallery holds no Session "${id}"`);
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
