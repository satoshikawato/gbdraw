// @ts-check
/** @import { LegendFontFaceTable } from '../utils/legend-font-metrics.generated.js' */
/** @import { RegisteredLegendFace } from './legend-layout.js' */
/** @import { FontFamilyList, FontFamilyName, GenericFontFamily } from './font-families.js' */
import { sha256Hex } from './byte-utils.js';
import { comparableFamilyName, parseFontFamilyList } from './font-families.js';
import { diagnosticError } from '../utils/error-normalization.js';

// The fonts registered in this page (custom-fonts D-01): the one record of
// each font's bytes and of what Python read from them. An entry is keyed by
// the SHA-256 of its bytes. Entries live for the page only; a Session stores
// the font family name, never the bytes (D-02). The list itself is app state
// (`state.loadedFonts`); this module holds no state of its own.

/** The SHA-256 of a font file's bytes, 64 lower-case hex digits. @typedef {string} FontIdentity */
/** Python's measurement name of a registered face (`user-<16 hex>-<style>`). @typedef {string} FontMeasurementName */

/**
 * One registered font face.
 * @typedef {object} WebFontEntry
 * @property {FontIdentity} identity
 * @property {FontFamilyName} family The typographic family name Python read.
 * @property {number} weight The OS/2 weight class.
 * @property {boolean} italic
 * @property {GenericFontFamily} generic
 * @property {FontFamilyList} fontFamily The stored list for the family, as Python wrote it (`font_family_list_for`).
 * @property {FontMeasurementName} measurementName
 * @property {'file' | 'local'} source A loaded font file (A) or an installed font (B, Local Font Access).
 * @property {string} label The file name, or the installed face's full name.
 * @property {ArrayBuffer} bytes Kept for the Worker replay and SVG embedding (D-05).
 * @property {LegendFontFaceTable} metrics Python's `font_metrics_table` of the face.
 */

/**
 * A font file to register: its bytes and where it came from.
 * @typedef {object} WebFontSource
 * @property {ArrayBuffer} bytes
 * @property {'file' | 'local'} source
 * @property {string} label
 */

/**
 * A font source with the identity of its bytes.
 * @typedef {WebFontSource & { identity: FontIdentity }} IdentifiedFontSource
 */

/**
 * What the Worker's `registerFonts` helper returns per font, in request order.
 * @typedef {object} RegisteredFontFacts
 * @property {FontIdentity} identity
 * @property {string} family
 * @property {FontFamilyList} fontFamily
 * @property {number} weight
 * @property {boolean} italic
 * @property {string} generic
 * @property {string} measurementName
 * @property {LegendFontFaceTable} metrics
 */

/**
 * The Worker's registration of font files (the `registerFonts` helper); its
 * reply is unvalidated (`RegisteredFontFacts[]` when well formed).
 * @typedef {(fonts: readonly { identity: FontIdentity, label: string, bytes: ArrayBuffer }[]) => Promise<unknown>} RegisterFontsInWorker
 */

/** @param {ArrayBuffer} bytes @returns {Promise<FontIdentity>} */
export const fontIdentityOf = (bytes) => sha256Hex(bytes);

/** @param {unknown} value @returns {value is RegisteredFontFacts} */
const isFontFacts = (value) => {
  if (!value || typeof value !== 'object') return false;
  const facts = /** @type {Record<string, unknown>} */ (value);
  return typeof facts.identity === 'string' && typeof facts.family === 'string' && facts.family.trim() !== ''
    && typeof facts.fontFamily === 'string' && facts.fontFamily.trim() !== ''
    && Number.isInteger(facts.weight) && typeof facts.italic === 'boolean'
    && ['serif', 'sans-serif', 'monospace'].includes(String(facts.generic))
    && typeof facts.measurementName === 'string' && Boolean(facts.metrics) && typeof facts.metrics === 'object';
};

// Local Font Access style names of the four faces a drawing uses (normal and
// bold weight, upright and italic), in the order they are registered.
const STANDARD_STYLES = Object.freeze([
  ['regular', 'normal', 'book', 'roman', 'plain'],
  ['bold'],
  ['italic', 'oblique'],
  ['bolditalic', 'boldoblique']
]);

/**
 * The faces of an installed family to register: its Regular, Bold, Italic,
 * and Bold Italic faces when it has them (Python's CSS matching picks among
 * the registered faces), else its first face.
 * @template {{ style: string }} T
 * @param {readonly T[]} faces
 * @returns {T[]}
 */
export const standardInstalledFaces = (faces) => {
  /** @param {T} face */
  const styleKey = (face) => face.style.toLowerCase().replace(/[\s_-]+/g, '');
  /** @type {T[]} */
  const chosen = [];
  STANDARD_STYLES.forEach((names) => {
    const face = faces.find((candidate) => names.includes(styleKey(candidate)));
    if (face) chosen.push(face);
  });
  return chosen.length > 0 ? chosen : faces.slice(0, 1);
};

/**
 * The sources with the identities of their bytes, each identity once, in
 * source order.
 * @param {readonly WebFontSource[]} sources
 * @returns {Promise<IdentifiedFontSource[]>}
 */
export const identifyWebFontSources = async (sources) => {
  /** @type {Map<FontIdentity, IdentifiedFontSource>} */
  const byIdentity = new Map();
  for (const source of sources) {
    const identity = await fontIdentityOf(source.bytes);
    if (!byIdentity.has(identity)) byIdentity.set(identity, { ...source, identity });
  }
  return [...byIdentity.values()];
};

/**
 * Register font files with Python and return their entries, in source order.
 * The Worker request carries the store's bytes: `postMessage` makes the one
 * copy the Worker needs.
 * @param {readonly IdentifiedFontSource[]} sources Each identity once.
 * @param {RegisterFontsInWorker} registerInWorker
 * @returns {Promise<WebFontEntry[]>}
 */
export const registerWebFonts = async (sources, registerInWorker) => {
  if (sources.length === 0) return [];
  const requested = sources.map(({ identity, label, bytes }) => ({ identity, label, bytes }));
  /** @type {unknown} */
  const reply = await registerInWorker(requested);
  /** @type {unknown[]} */
  const facts = Array.isArray(reply) ? reply : [];
  if (facts.length !== requested.length) {
    throw diagnosticError('RESULT_INVALID', { field: 'font_family' }, { stage: 'result-admission', operation: 'registerFonts' });
  }
  return facts.map((fact, index) => {
    if (!isFontFacts(fact) || fact.identity !== requested[index].identity) {
      throw diagnosticError('RESULT_INVALID', { field: 'font_family' }, { stage: 'result-admission', operation: 'registerFonts' });
    }
    const source = sources[index];
    return Object.freeze({
      identity: fact.identity,
      family: fact.family,
      weight: fact.weight,
      italic: fact.italic,
      generic: fact.generic,
      fontFamily: fact.fontFamily,
      measurementName: fact.measurementName,
      source: source.source,
      label: source.label,
      bytes: source.bytes,
      metrics: fact.metrics
    });
  });
};

/**
 * The browser's check of each font before Python registers it: a FontFace
 * per identity, loaded under a placeholder family until Python names it.
 * @param {readonly IdentifiedFontSource[]} sources
 * @param {typeof FontFace} FontFaceClass
 * @returns {Promise<Map<FontIdentity, FontFace>>}
 */
const loadDocumentFaces = async (sources, FontFaceClass) => {
  /** @type {Map<FontIdentity, FontFace>} */
  const faces = new Map();
  for (const source of sources) {
    // FontFace copies the bytes it is given.
    const face = new FontFaceClass(`gbdraw-font-${source.identity}`, source.bytes);
    try {
      await face.load();
    } catch (_error) {
      throw diagnosticError('INPUT_INVALID', { field: 'font_family', reason: 'FONT_FILE' }, { stage: 'font-validation', operation: 'registerFonts' });
    }
    faces.set(source.identity, face);
  }
  return faces;
};

/**
 * The result of loading font sources into the page.
 * @typedef {object} LoadedWebFonts
 * @property {WebFontEntry[]} added The entries registered by this load.
 * @property {WebFontEntry[]} chosen The entry of each source's identity, in source order (stored or added).
 */

/**
 * Load font sources into the page: the browser checks each new font first
 * (FontFace), then Python registers it, then the preview gets its face, so
 * the store, Python, and the preview never disagree (OV-687). A font already
 * in the store (by the SHA-256 of its bytes) is not registered or added again.
 * @param {readonly WebFontSource[]} sources
 * @param {object} options
 * @param {readonly WebFontEntry[]} options.known The page's fonts.
 * @param {RegisterFontsInWorker} options.registerInWorker
 * @param {{ add: (face: FontFace) => unknown } | null} [options.fontSet] The document's fonts.
 * @param {typeof FontFace | null} [options.FontFaceClass]
 * @returns {Promise<LoadedWebFonts>}
 */
export const loadWebFonts = async (sources, {
  known,
  registerInWorker,
  fontSet = globalThis.document?.fonts ?? null,
  FontFaceClass = typeof FontFace === 'function' ? FontFace : null
}) => {
  const identified = await identifyWebFontSources(sources);
  const stored = new Map(known.map((entry) => [entry.identity, entry]));
  const fresh = identified.filter((source) => !stored.has(source.identity));
  const faces = fontSet && FontFaceClass ? await loadDocumentFaces(fresh, FontFaceClass) : null;
  const added = await registerWebFonts(fresh, registerInWorker);
  added.forEach((entry) => {
    const face = faces?.get(entry.identity);
    if (!face || !fontSet) return;
    face.family = entry.family;
    face.weight = String(entry.weight);
    face.style = entry.italic ? 'italic' : 'normal';
    fontSet.add(face);
  });
  added.forEach((entry) => stored.set(entry.identity, entry));
  return {
    added,
    chosen: identified.map((source) => /** @type {WebFontEntry} */ (stored.get(source.identity)))
  };
};

/**
 * The page's fonts after `added` were registered: an added face replaces an
 * entry with its identity; the rest keep their order.
 * @param {readonly WebFontEntry[]} current
 * @param {readonly WebFontEntry[]} added
 * @returns {WebFontEntry[]}
 */
export const mergeWebFonts = (current, added) => {
  const addedIdentities = new Set(added.map((entry) => entry.identity));
  return [...current.filter((entry) => !addedIdentities.has(entry.identity)), ...added];
};

/**
 * The families of the page's fonts, each once, in registration order.
 * @param {readonly WebFontEntry[]} entries
 * @returns {{ family: FontFamilyName, fontFamily: FontFamilyList, source: 'file' | 'local' }[]}
 */
export const webFontFamilies = (entries) => {
  /** @type {Map<string, { family: FontFamilyName, fontFamily: FontFamilyList, source: 'file' | 'local' }>} */
  const families = new Map();
  entries.forEach((entry) => {
    const key = comparableFamilyName(entry.family);
    if (!families.has(key)) families.set(key, { family: entry.family, fontFamily: entry.fontFamily, source: entry.source });
  });
  return [...families.values()];
};

/**
 * The names of the page's font families, each once (the registration state
 * the custom-font and availability rules of font-families.js take).
 * @param {readonly WebFontEntry[]} entries
 * @returns {FontFamilyName[]}
 */
export const webFontFamilyNames = (entries) => webFontFamilies(entries).map((entry) => entry.family);

/**
 * The faces the Legend layout measures with (`withRegisteredFontFaces`).
 * @param {readonly WebFontEntry[]} entries
 * @returns {RegisteredLegendFace[]}
 */
export const registeredLegendFaces = (entries) => entries.map((entry) => ({
  family: entry.family,
  weight: entry.weight,
  italic: entry.italic,
  measurementName: entry.measurementName,
  table: entry.metrics
}));

/**
 * The identities of the page's fonts a `font-family` list names: what a
 * Generate needs registered in the Worker.
 * @param {readonly WebFontEntry[]} entries
 * @param {unknown} fontFamily
 * @returns {FontIdentity[]}
 */
export const webFontIdentitiesFor = (entries, fontFamily) => {
  const named = new Set(parseFontFamilyList(fontFamily).map(comparableFamilyName));
  return entries.filter((entry) => named.has(comparableFamilyName(entry.family))).map((entry) => entry.identity);
};

/**
 * The text font a canonical render request names (`objects.text.font_family`,
 * as an override or in its config), or null: what a Generate needs registered
 * and what the shown Result was drawn with.
 * @param {unknown} request
 * @returns {unknown}
 */
export const requestFontFamily = (request) => {
  const options = /** @type {{ diagramOptions?: { configOverrides?: Record<string, unknown>, config?: { objects?: { text?: { font_family?: unknown } } } } } | null} */ (request)?.diagramOptions;
  return options?.configOverrides?.['objects.text.font_family'] ?? options?.config?.objects?.text?.font_family ?? null;
};
