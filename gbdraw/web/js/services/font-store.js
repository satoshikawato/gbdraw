// @ts-check
/** @import { LegendFontFaceTable } from '../utils/legend-font-metrics.generated.js' */
/** @import { RegisteredLegendFace } from './legend-layout.js' */
/** @import { FontFamilyList, FontFamilyName, GenericFontFamily } from './font-families.js' */
import { comparableFamilyName, fontFamilyListFor, parseFontFamilyList } from './font-families.js';
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
 * @property {FontFamilyList} fontFamily The stored list for the family (`font_family_list_for`).
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
 * What the Worker's `registerFonts` helper returns per font, in request order.
 * @typedef {object} RegisteredFontFacts
 * @property {FontIdentity} identity
 * @property {string} family
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
export const fontIdentityOf = async (bytes) => {
  const digest = new Uint8Array(await globalThis.crypto.subtle.digest('SHA-256', bytes));
  return Array.from(digest, (value) => value.toString(16).padStart(2, '0')).join('');
};

/** @param {unknown} value @returns {value is RegisteredFontFacts} */
const isFontFacts = (value) => {
  if (!value || typeof value !== 'object') return false;
  const facts = /** @type {Record<string, unknown>} */ (value);
  return typeof facts.identity === 'string' && typeof facts.family === 'string' && facts.family.trim() !== ''
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
 * Register font files with Python and return their entries, in source order.
 * Identical bytes are registered once.
 * @param {readonly WebFontSource[]} sources
 * @param {RegisterFontsInWorker} registerInWorker
 * @returns {Promise<WebFontEntry[]>}
 */
export const registerWebFonts = async (sources, registerInWorker) => {
  /** @type {Map<FontIdentity, WebFontSource>} */
  const byIdentity = new Map();
  for (const source of sources) {
    const identity = await fontIdentityOf(source.bytes);
    if (!byIdentity.has(identity)) byIdentity.set(identity, source);
  }
  const requested = Array.from(byIdentity, ([identity, source]) => ({ identity, label: source.label, bytes: source.bytes.slice(0) }));
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
    const source = /** @type {WebFontSource} */ (byIdentity.get(fact.identity));
    return Object.freeze({
      identity: fact.identity,
      family: fact.family,
      weight: fact.weight,
      italic: fact.italic,
      generic: fact.generic,
      fontFamily: fontFamilyListFor(fact.family, fact.generic),
      measurementName: fact.measurementName,
      source: source.source,
      label: source.label,
      bytes: source.bytes,
      metrics: fact.metrics
    });
  });
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
 * Make each entry drawable in the inline preview, which uses the document's fonts.
 * @param {readonly WebFontEntry[]} entries
 * @param {FontFaceSet | null | undefined} [fontSet]
 * @returns {Promise<void>}
 */
export const addWebFontsToDocument = async (entries, fontSet = globalThis.document?.fonts) => {
  if (!fontSet || typeof FontFace !== 'function') return;
  for (const entry of entries) {
    const face = new FontFace(entry.family, entry.bytes.slice(0), {
      weight: String(entry.weight),
      style: entry.italic ? 'italic' : 'normal'
    });
    await face.load();
    fontSet.add(face);
  }
};
