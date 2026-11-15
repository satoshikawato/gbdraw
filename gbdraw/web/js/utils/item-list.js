// @ts-check
// Shown items one after the other with no conjunction (record ids, gene names, file names, ranges, translated
// names), separated as the shown UI language separates them: the one owner of that separator. A list that ends in
// "and" is `listText` in i18n.js.
/** @import { UiLocale } from './i18n.js' */
import { uiLocale } from './i18n.js';
import { ownEntry } from './own-entry.js';

// Right-to-Left Isolate, First Strong Isolate and Pop Directional Isolate (Unicode bidi).
const RLI = '⁧';
const FSI = '⁨';
const PDI = '⁩';

/**
 * The item separator of each UI language that reads right to left and writes the Arabic comma ، (U+060C) between
 * list items; every other language writes ", " (I18N_STYLE: Persian and Arabic punctuation).
 * @type {Readonly<Partial<Record<UiLocale, string>>>}
 */
const RIGHT_TO_LEFT_ITEM_SEPARATORS = Object.freeze({ ar: '، ', fa: '، ' });

/**
 * Shown items as one list: `A, B, C` in English and every left-to-right language (the same bytes as
 * `items.join(', ')`), `A، B، C` in Persian and Arabic. There each item is isolated (FSI ... PDI), so a record id or a
 * gene name keeps its own order, and the list reads right to left as one isolate (RLI ... PDI) wherever the text puts
 * it. No items show as the empty text. Reading it makes a render depend on the UI language.
 * @param {readonly string[]} items each already a shown text
 * @returns {string}
 */
export const itemListText = (items) => {
  const separator = ownEntry(RIGHT_TO_LEFT_ITEM_SEPARATORS, uiLocale.value);
  if (separator === undefined || items.length === 0) return items.join(', ');
  return `${RLI}${items.map((item) => `${FSI}${item}${PDI}`).join(separator)}${PDI}`;
};
