// The bytes a Session Save writes for each input, as base64 descriptors, so a
// test reads what the Session file holds: a file-backed descriptor (a File
// Generate described, or an entry a Load copied) is read through
// `storedResourceBytes`, the one owner of what Save stores (a gzip input as it
// is, with its size); any other descriptor is written as it holds its bytes.
// `backing` is the session-resource-backing module the test's descriptors were made by.
export const writtenResourcesWith = ({ fileBackedReaderOf, storedResourceBytes }) => async (resources) => Object.fromEntries(
  await Promise.all(Object.entries(resources).map(async ([id, descriptor]) => {
    if (!fileBackedReaderOf(descriptor)) return [id, descriptor];
    const stored = storedResourceBytes(descriptor);
    const bytes = Buffer.from(await new Response(await stored.open()).arrayBuffer());
    return [id, { ...descriptor, size: stored.size, encoding: 'base64', data: bytes.toString('base64') }];
  }))
);

/**
 * A Gallery Session as a Session file container, its document changed by
 * `edit` first.
 * @param {string} id
 * @param {(document: Record<string, any>) => void} [edit]
 * @returns {Promise<Blob>}
 */
export const gallerySessionFile = async (id, edit = () => {}) => {
  const { readFileSync } = await import('node:fs');
  const sessions = new URL('../../../gbdraw/web/gallery/sessions/', import.meta.url);
  const { document, bytes } = await readSessionFile(new Blob([readFileSync(new URL(`${id}.gbdraw`, sessions))]));
  edit(document);
  return writeSessionFile(document, bytes);
};

/**
 * A Session file container's document and the bytes of each resource entry.
 * @param {Blob} blob
 * @returns {Promise<{ document: Record<string, any>, bytes: Record<string, Uint8Array> }>}
 */
export const readSessionFile = async (blob) => {
  const { openSessionContainerEntry, readSessionContainer } = await import('../../../gbdraw/web/js/services/session-container.js');
  const read = await readSessionContainer(blob);
  /** @type {Record<string, Uint8Array>} */
  const bytes = {};
  for (const [resourceId, entry] of Object.entries(read.resources)) {
    bytes[resourceId] = new Uint8Array(await new Response(openSessionContainerEntry(blob, entry)).arrayBuffer());
  }
  return { document: /** @type {Record<string, any>} */ (read.document), bytes };
};

/**
 * A Session file container of `document` whose resource entries hold `bytes`.
 * @param {Record<string, any>} document
 * @param {Record<string, Uint8Array>} bytes
 * @returns {Promise<Blob>}
 */
export const writeSessionFile = async (document, bytes) => {
  const { writeSessionContainer } = await import('../../../gbdraw/web/js/services/session-container.js');
  return writeSessionContainer({
    document,
    documentBytes: () => new Blob([JSON.stringify(document)]).stream(),
    readResource: (resourceId) => new Blob([bytes[resourceId]]).stream()
  });
};
