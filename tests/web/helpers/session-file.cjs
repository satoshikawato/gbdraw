'use strict';

// Session files as specs read and write them: a saved Session is a Session
// file container (`.gbdraw`), read here with the Web app's own container
// reader, never by unzipping it another way.
const fs = require('node:fs/promises');
const path = require('node:path');
const { pathToFileURL } = require('node:url');

const containerModule = () => import(pathToFileURL(
  path.join(__dirname, '../../../gbdraw/web/js/services/session-container.js')).href);
const writtenResourcesModule = () => import(pathToFileURL(path.join(__dirname, 'written-resources.mjs')).href);
const sessionDocumentsModule = () => import(pathToFileURL(path.join(__dirname, 'session-documents.mjs')).href);

/** The Session document a saved `.gbdraw` file holds. */
const readSessionFile = async (file) => {
  const { readSessionContainer } = await containerModule();
  return (await readSessionContainer(new Blob([await fs.readFile(file)]))).document;
};

/**
 * A saved `.gbdraw` file's document and, for each resource it declares, its
 * entry's compression method (0 stored, 8 deflated) and bytes.
 */
const readSessionFileEntries = async (file) => {
  const { openSessionContainerEntry, readSessionContainer } = await containerModule();
  const blob = new Blob([await fs.readFile(file)]);
  const container = await readSessionContainer(blob);
  const resources = {};
  for (const [id, entry] of Object.entries(container.resources)) {
    resources[id] = { method: entry.method, bytes: Buffer.from(await new Response(openSessionContainerEntry(blob, entry)).arrayBuffer()) };
  }
  return { document: container.document, resources };
};

/** The bytes of the Gallery Session `id` as a Session file container, its document changed by `edit` first. */
const gallerySessionBytes = async (id, edit = () => {}) => {
  const { gallerySessionFile } = await writtenResourcesModule();
  return Buffer.from(await (await gallerySessionFile(id, edit)).arrayBuffer());
};

/** Write the Gallery Session `id` as a Session file container at `file`; returns `file`. */
const writeGallerySessionFile = async (id, file, edit = () => {}) => {
  await fs.writeFile(file, await gallerySessionBytes(id, edit));
  return file;
};

/**
 * A saved Session file (a `.gbdraw` file, or an older JSON Session) in the
 * JSON form a test edits: each container entry back as a base64 descriptor
 * (`readSavedSession`).
 */
const readSessionFileJson = async (file) => (await sessionDocumentsModule()).readSavedSession(new Blob([await fs.readFile(file)]));

/** The bytes of the Session file of a JSON-form `document` (`sessionFileOf`): a `.gbdraw` file for a Session 51 or later. */
const sessionFileBytesOf = async (document) => Buffer.from(await (await (await sessionDocumentsModule()).sessionFileOf(document)).arrayBuffer());

module.exports = {
  gallerySessionBytes, readSessionFile, readSessionFileEntries, readSessionFileJson, sessionFileBytesOf, writeGallerySessionFile
};
