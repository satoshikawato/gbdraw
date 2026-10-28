// @ts-check
// The Interactive Gallery's addresses and tutorial state (#1383, PD-OI-108).
// The Gallery page (`gallery/gallery.js`) and the app's tutorial links read
// them here, so an example's page and whether its tutorial is published have
// one owner.

/**
 * The Gallery that the header and the footer open. The pip GUI ships no
 * Gallery, so every build opens gbdraw.app's.
 */
export const INTERACTIVE_GALLERY_URL = 'https://gbdraw.app/gallery/';

/** @typedef {'preview' | 'tutorial' | 'command' | 'files'} GalleryTab */

/**
 * The hash that opens each tab of an example's page.
 * @type {Readonly<Record<GalleryTab, string>>}
 */
export const GALLERY_TAB_HASHES = Object.freeze({
  preview: '',
  tutorial: 'Tutorial',
  command: 'Command',
  files: 'Files'
});

/**
 * The path of an example's page: `<Gallery root path>/<id>`.
 * @param {string} galleryRootPath the Gallery's path without a final slash, e.g. `/gallery`
 * @param {string} exampleId
 * @returns {string}
 */
export const galleryExamplePath = (galleryRootPath, exampleId) => `${galleryRootPath}/${encodeURIComponent(exampleId)}`;

/** @typedef {'ready' | 'draft' | 'planned'} GalleryTutorialStatus */
/** @typedef {{ tutorial?: unknown, tutorialStatus?: unknown }} GalleryTutorialEntry */

/**
 * The state of an example's tutorial in the catalog (`examples.json`): its
 * `tutorialStatus`, or `ready` for an entry that names tutorial data without one.
 * @param {GalleryTutorialEntry | null | undefined} entry
 * @returns {GalleryTutorialStatus}
 */
export const galleryTutorialStatus = (entry) => {
  const status = typeof entry?.tutorialStatus === 'string' ? entry.tutorialStatus.trim().toLowerCase() : '';
  if (status === 'ready' || status === 'draft' || status === 'planned') return status;
  return entry?.tutorial ? 'ready' : 'planned';
};

/**
 * Whether the public Gallery shows the example's tutorial: the entry names its
 * tutorial data and the tutorial is ready.
 * @param {GalleryTutorialEntry | null | undefined} entry
 * @returns {boolean}
 */
export const galleryTutorialPublished = (entry) => Boolean(entry?.tutorial) && galleryTutorialStatus(entry) === 'ready';
