// @ts-check
/** @import { FontFamilyName } from './font-families.js' */
import { t, tn } from '../utils/i18n.js';
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
const nameList = (names) => {
  if (names.length < 2) return names.join('');
  if (names.length === 2) return t('{first} and {second}', { first: names[0], second: names[1] });
  return t('{list}, and {last}', { list: names.slice(0, -1).join(', '), last: String(names.at(-1)) });
};

/** @template T @param {readonly T[]} values @returns {T[]} */
const unique = (values) => [...new Set(values)];

/**
 * The Save Session notice's warnings: one per custom family of the drawings'
 * fonts (`customFontFamilies`).
 * @param {readonly unknown[]} fontFamilies Each drawing's `adv.font_family`.
 * @param {Iterable<FontFamilyName>} registeredFamilies The families loaded in this page.
 * @returns {string[]}
 */
export const customFontSaveWarnings = (fontFamilies, registeredFamilies) => {
  const registered = [...registeredFamilies];
  return unique(fontFamilies.flatMap((fontFamily) => customFontFamilies(fontFamily, registered)))
    .map((family) => t('The Session keeps the font name {family}, not the font file. To open it elsewhere, the same font must be available there, or it is replaced with the nearest available font.', { family }));
};

/**
 * The fonts a Session names that this page cannot draw where the family drawn
 * differs from the requested first family (Python's `replaced`, OV-743), and
 * the nearest font each such list is drawn with instead (the first family of
 * the list Python resolves to); null when every list draws its first family.
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
    .filter((resolution) => resolution.replaced.length > 0);
  if (resolutions.length === 0) return null;
  return {
    missing: unique(resolutions.flatMap((resolution) => resolution.replaced)),
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
  const names = nameList(missing);
  const substitutes = nameList(nearest);
  return {
    title: t('This Session uses a font that is not available'),
    description: tn(missing.length,
      'Its text is set in {names}, which is not part of gbdraw and is not loaded in this page.',
      'Its text is set in {names}, which are not part of gbdraw and are not loaded in this page.', { names }),
    items: loadedFamily ? [t('The loaded font is {loaded}, not {names}.', { loaded: loadedFamily, names })] : [],
    note: tn(nearest.length,
      'Load with {substitutes}: the Session opens and its text is drawn in {substitutes}, the nearest available font.',
      'Load with {substitutes}: the Session opens and its text is drawn in {substitutes}, the nearest available fonts.', { substitutes })
      + ' ' + tn(missing.length,
        'The Font setting keeps the name {names}, so the font is used again once it is loaded.',
        'The Font setting keeps the names {names}, so the fonts are used again once they are loaded.', { names })
      + ' ' + t('The saved preview keeps its look until you Generate.'),
    actions: [
      { id: SESSION_FONT_ACTIONS.FILE, label: t('Load font file…'), tone: /** @type {const} */ ('primary') },
      ...(installedFontsSupported
        ? [{ id: SESSION_FONT_ACTIONS.INSTALLED, label: t('Use installed font…'), tone: /** @type {const} */ ('secondary') }]
        : []),
      { id: SESSION_FONT_ACTIONS.SUBSTITUTE, label: t('Load with {substitutes}', { substitutes }), tone: /** @type {const} */ ('secondary') }
    ],
    cancelLabel: t('Cancel')
  };
};
