// The bytes a Session Save writes for each input, as base64 descriptors, so a
// test reads what the Session file holds: a file-backed descriptor (a File
// Generate described, or an entry a Load copied) is read through
// `sessionResourceStream`, the one reader Save streams an entry from; any other
// descriptor is written as it holds its bytes. `backing` is the
// session-resource-backing module the test's descriptors were made by.
export const writtenResourcesWith = ({ fileBackedReaderOf, sessionResourceStream }) => async (resources) => Object.fromEntries(
  await Promise.all(Object.entries(resources).map(async ([id, descriptor]) => {
    if (!fileBackedReaderOf(descriptor)) return [id, descriptor];
    const bytes = Buffer.from(await new Response(await sessionResourceStream(descriptor)).arrayBuffer());
    return [id, { ...descriptor, encoding: 'base64', data: bytes.toString('base64') }];
  }))
);

/**
 * A Gallery Session as a Session file container, its document changed by
 * `edit` first: the Gallery's `<id>.gbdraw` when it holds one, else its JSON
 * Session with each base64 input written as an entry.
 * @param {string} id
 * @param {(document: Record<string, any>) => void} [edit]
 * @returns {Promise<Blob>}
 */
export const gallerySessionFile = async (id, edit = () => {}) => {
  const { existsSync, readFileSync } = await import('node:fs');
  const { openSessionContainerEntry, readSessionContainer, writeSessionContainer } = await import(
    '../../../gbdraw/web/js/services/session-container.js');
  const sessions = new URL('../../../gbdraw/web/gallery/sessions/', import.meta.url);
  const container = new URL(`${id}.gbdraw`, sessions);
  /** @type {Record<string, any>} */
  let document;
  /** @type {Record<string, Uint8Array>} */
  const bytes = {};
  if (existsSync(container)) {
    const blob = new Blob([readFileSync(container)]);
    const read = await readSessionContainer(blob);
    document = /** @type {Record<string, any>} */ (read.document);
    for (const [resourceId, entry] of Object.entries(read.resources)) {
      bytes[resourceId] = new Uint8Array(await new Response(openSessionContainerEntry(blob, entry)).arrayBuffer());
    }
  } else {
    document = JSON.parse(readFileSync(new URL(`${id}.gbdraw-session.json`, sessions), 'utf8'));
    for (const [resourceId, descriptor] of Object.entries(document.resources)) {
      const { encoding: _encoding, data, checksum: _checksum, ...fields } = descriptor;
      bytes[resourceId] = Buffer.from(data, 'base64');
      document.resources[resourceId] = { kind: 'web-file', type: 'application/octet-stream', lastModified: 0, ...fields };
    }
  }
  edit(document);
  return writeSessionContainer({
    document,
    documentBytes: () => new Blob([JSON.stringify(document)]).stream(),
    readResource: (resourceId) => new Blob([bytes[resourceId]]).stream()
  });
};
