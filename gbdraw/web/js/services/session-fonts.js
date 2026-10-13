// @ts-check
/** @import { FontFamilyName } from './font-families.js' */
import { customFontFamilies, isGenericFamily, parseFontFamilyList, resolveWebFontFamily } from './font-families.js';

// What a Session says about its fonts (custom fonts D-02, D-03; Product
// Contract diagram-generation.diagram-font rev 1, 〔警告〕 and 〔読み込むとき〕):
// Save warns for each font a Session keeps by name only, and Load asks before
// a Session whose font is neither bundled nor loaded in this page replaces the
// current one. Both decide with the Python-pinned rules of font-families.js,
// in the browser, so the check constructs no diagram Worker.

/** The ids of the Load question's choices. */
export const SESSION_FONT_ACTIONS = Object.freeze({
  FILE: 'file',
  INSTALLED: 'installed',
  SUBSTITUTE: 'substitute'
});

/** @param {readonly string[]} names "A", "A and B", "A, B, and C" */
const nameList = (names) => (names.length < 3
  ? names.join(' and ')
  : `${names.slice(0, -1).join(', ')}, and ${names.at(-1)}`);

/** @template T @param {readonly T[]} values @returns {T[]} */
const unique = (values) => [...new Set(values)];

/**
 * The Save Session notice's warnings: one per family of the drawings' fonts
 * that gbdraw does not bundle (`customFontFamilies`).
 * @param {readonly unknown[]} fontFamilies Each drawing's `adv.font_family`.
 * @returns {string[]}
 */
export const customFontSaveWarnings = (fontFamilies) => unique(fontFamilies.flatMap(customFontFamilies))
  .map((family) => `The Session keeps the font name ${family}, not the font file. To open it elsewhere, the same font must be available there, or it is replaced with the nearest available font.`);

/**
 * The fonts a Session names that this page cannot draw, and the nearest font
 * each list is drawn with instead (the first family of the list Python
 * resolves to); null when every family resolves.
 * @typedef {{ missing: FontFamilyName[], nearest: FontFamilyName[] }} MissingSessionFonts
 */

/**
 * @param {readonly unknown[]} fontFamilies Each drawing's `adv.font_family` in the Session.
 * @param {Iterable<FontFamilyName>} registeredFamilies The families loaded in this page.
 * @returns {MissingSessionFonts | null}
 */
export const missingSessionFonts = (fontFamilies, registeredFamilies) => {
  const registered = [...registeredFamilies];
  const resolutions = fontFamilies
    .filter((value) => typeof value === 'string' && value.trim() !== '')
    .map((value) => resolveWebFontFamily(String(value), registered))
    .filter((resolution) => resolution.missing.length > 0);
  if (resolutions.length === 0) return null;
  return {
    missing: unique(resolutions.flatMap((resolution) => resolution.missing)),
    nearest: unique(resolutions.map((resolution) => (
      parseFontFamilyList(resolution.effective).find((family) => !isGenericFamily(family)) || resolution.effective
    )))
  };
};

/**
 * The Load question's request for the app's choice dialog.
 * @param {MissingSessionFonts} fonts
 * @param {{ installedFontsSupported: boolean, loadedFamily: FontFamilyName | null }} options
 *   `loadedFamily`: the family of the font just loaded, when it is not the one the Session names.
 */
export const sessionFontQuestion = ({ missing, nearest }, { installedFontsSupported, loadedFamily }) => {
  const one = missing.length === 1;
  const names = nameList(missing);
  const substitutes = nameList(nearest);
  return {
    title: 'This Session uses a font that is not available',
    description: `Its text is set in ${names}, which ${one ? 'is' : 'are'} not part of gbdraw and ${one ? 'is' : 'are'} not loaded in this page.`,
    items: loadedFamily ? [`The loaded font is ${loadedFamily}, not ${names}.`] : [],
    note: `Load with ${substitutes}: the Session opens and its text is drawn in ${substitutes}, `
      + `the nearest available font${nearest.length === 1 ? '' : 's'}. `
      + (one
        ? `The Font setting keeps the name ${names}, so the font is used again once it is loaded. `
        : `The Font setting keeps the names ${names}, so the fonts are used again once they are loaded. `)
      + 'The saved preview keeps its look until you Generate.',
    actions: [
      { id: SESSION_FONT_ACTIONS.FILE, label: 'Load font file…', tone: /** @type {const} */ ('primary') },
      ...(installedFontsSupported
        ? [{ id: SESSION_FONT_ACTIONS.INSTALLED, label: 'Use installed font…', tone: /** @type {const} */ ('secondary') }]
        : []),
      { id: SESSION_FONT_ACTIONS.SUBSTITUTE, label: `Load with ${substitutes}`, tone: /** @type {const} */ ('secondary') }
    ],
    cancelLabel: 'Cancel'
  };
};
