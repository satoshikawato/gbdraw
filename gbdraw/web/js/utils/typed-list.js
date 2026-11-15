// @ts-check

// The separators of a list typed into one text field: ring names and the Label filtering
// Blacklist (its text and the words of a Blacklist file). A comma, the Arabic comma ،
// (U+060C), which fa and ar keyboards type and their catalogs write for "comma", or a line
// break (CRLF, LF or CR) separates two items (OV-2792). No typed list takes ';', so the
// Arabic semicolon is no separator either. Text the user does not type as a list (slot specs,
// numbers, colors, font stacks, joins the Web writes) keeps its own ',' rules; Python reads
// these lists as arrays the Web sends and never splits this text.

/** One separator: CRLF is one line break. */
const SEPARATOR = /(\r\n|[\r\n,،])/;
/** A run of separators between two items. */
const SEPARATOR_RUN = /[\r\n,،]+/;
/** A separator at the end of a text. */
const TRAILING_SEPARATOR = /[\r\n,،]$/;

/**
 * The items of a typed list: the text cut at each run of separators, each item trimmed, blank
 * items dropped. A missing or empty value has none.
 * @param {unknown} text
 * @returns {string[]}
 */
export const splitTypedList = (text) => String(text || '')
  .split(SEPARATOR_RUN)
  .map((item) => item.trim())
  .filter(Boolean);

/**
 * A typed list's text cut at each separator, with the separators kept: the items (untrimmed,
 * blank ones too) at even places and the one separator after each at the odd places, so
 * `join('')` gives the text back.
 * @param {string} text
 * @returns {string[]}
 */
export const splitTypedListKeepingSeparators = (text) => text.split(SEPARATOR);

/**
 * Whether a typed list's text ends with a separator, so an item appended needs none.
 * @param {string} text
 */
export const endsWithTypedListSeparator = (text) => TRAILING_SEPARATOR.test(text);
