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

/** The Session document a saved `.gbdraw` file holds. */
const readSessionFile = async (file) => {
  const { readSessionContainer } = await containerModule();
  return (await readSessionContainer(new Blob([await fs.readFile(file)]))).document;
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

module.exports = { gallerySessionBytes, readSessionFile, writeGallerySessionFile };
