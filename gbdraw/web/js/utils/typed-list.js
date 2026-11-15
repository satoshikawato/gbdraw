// @ts-check

// The separators of a list typed into one text field: ring names and the Label filtering
// Blacklist (its text and the words of a Blacklist file), OV-2792. Python owns them
// (`TYPED_LIST_SEPARATORS` in gbdraw/core/typed_list.py, which its CLI options and stored
// `label_blacklist` text read too), and the contract generator carries them here: the comma,
// the Arabic comma ، (U+060C) that fa and ar keyboards type and their catalogs write for
// "comma", and a line break, CRLF as one. Text the user does not type as a list (slot specs,
// numbers, colors, font stacks, joins the Web writes) keeps its own ',' rules.
import { TYPED_LIST_SEPARATORS } from '../contracts/schema.generated.js';

/** The separators as alternatives of a pattern, CRLF before CR and LF as Python orders them. */
const SEPARATORS = TYPED_LIST_SEPARATORS.map((separator) => separator.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')).join('|');
/** One separator. */
const SEPARATOR = new RegExp(`(${SEPARATORS})`);
/** A run of separators between two items. */
const SEPARATOR_RUN = new RegExp(`(?:${SEPARATORS})+`);
/** A separator at the end of a text. */
const TRAILING_SEPARATOR = new RegExp(`(?:${SEPARATORS})$`);

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
