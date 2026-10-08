// @ts-check
// UJ-06 (Owner 2026-10-08): the Gallery catalog that this origin serves is the
// one source of the empty state's examples. gbdraw.app serves it; the pip GUI
// does not ship it, so there the list is empty and Load an example is absent.

/** @typedef {{ text: string, italic: boolean }} GalleryTitlePart */
/**
 * @typedef {{
 *   id: string,
 *   title: string,
 *   titleParts: GalleryTitlePart[],
 *   mode: 'Circular' | 'Linear' | '',
 *   sessionUrl: string,
 *   sessionName: string,
 *   thumbnailUrl: string
 * }} GalleryExample
 */

const GALLERY_CATALOG_URL = new URL('../../gallery/examples.json', import.meta.url);
const SESSION_FILE_NAME = /\.json(?:\.gz)?$/i;
const MODES = /** @type {const} */ (['Circular', 'Linear']);

/**
 * A catalog title marks taxon names with `<i>`; every other character is text.
 * @param {string} title
 * @returns {GalleryTitlePart[]}
 */
export const galleryTitleParts = (title) => {
  /** @type {GalleryTitlePart[]} */
  const parts = [];
  let italic = false;
  for (const token of title.split(/(<\/?i>)/)) {
    if (token === '<i>' || token === '</i>') italic = token === '<i>';
    else if (token) parts.push({ text: token, italic });
  }
  return parts;
};

/**
 * A same-origin URL that `reference` names relative to the catalog, or null.
 * @param {unknown} reference
 * @param {URL} catalogUrl
 * @returns {URL | null}
 */
const sameOriginUrl = (reference, catalogUrl) => {
  if (typeof reference !== 'string' || !reference) return null;
  try {
    const url = new URL(reference, catalogUrl);
    return url.origin === catalogUrl.origin ? url : null;
  } catch {
    return null;
  }
};

/**
 * The usable entries of a catalog, in catalog order. An entry is usable when it
 * names an id, a title, and a same-origin `.json` or `.json.gz` Session.
 * @param {unknown} catalog
 * @param {URL} [catalogUrl]
 * @returns {GalleryExample[]}
 */
export const galleryExamplesFromCatalog = (catalog, catalogUrl = GALLERY_CATALOG_URL) => {
  if (!Array.isArray(catalog)) return [];
  /** @type {GalleryExample[]} */
  const examples = [];
  for (const entry of catalog) {
    if (!entry || typeof entry !== 'object') continue;
    const { id, title, session, thumbnail, tags, mode } = /** @type {Record<string, unknown>} */ (entry);
    const sessionUrl = sameOriginUrl(session, catalogUrl);
    const sessionName = decodeURIComponent(sessionUrl?.pathname.split('/').pop() || '');
    if (typeof id !== 'string' || !id || typeof title !== 'string' || !title.trim()
      || !sessionUrl || !SESSION_FILE_NAME.test(sessionName)) continue;
    const labels = [mode, ...(Array.isArray(tags) ? tags : [])]
      .map((value) => (typeof value === 'string' ? value.toLowerCase() : ''));
    examples.push(Object.freeze({
      id,
      title,
      titleParts: galleryTitleParts(title),
      mode: MODES.find((name) => labels.includes(name.toLowerCase())) || '',
      sessionUrl: sessionUrl.href,
      sessionName,
      thumbnailUrl: sameOriginUrl(thumbnail, catalogUrl)?.href || ''
    }));
  }
  return examples;
};

/**
 * Reads the catalog once. A missing or unreadable catalog gives no examples and
 * no error: the examples are optional.
 * @returns {Promise<GalleryExample[]>}
 */
export const readGalleryExamples = async () => {
  try {
    const response = await fetch(GALLERY_CATALOG_URL);
    return response.ok ? galleryExamplesFromCatalog(await response.json()) : [];
  } catch {
    return [];
  }
};

/**
 * The example's Session as a File named like its URL, so `.gz` is read as gzip.
 * @param {GalleryExample} example
 * @returns {Promise<File>}
 */
export const fetchGalleryExampleSession = async (example) => {
  const response = await fetch(example.sessionUrl);
  if (!response.ok) throw new Error(`The example Session could not be read (HTTP ${response.status}).`);
  return new File([await response.blob()], example.sessionName, {
    type: example.sessionName.toLowerCase().endsWith('.gz') ? 'application/gzip' : 'application/json'
  });
};
