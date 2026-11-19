// @ts-check
// The one display of user-facing error text (D-17). The wording owner,
// `normalizeUserFacingError` (R6), returns the English `summary` and the parts
// it is composed of (`summaryParts`; `noteParts` of a live-edit failure). A
// display shows the parts through `errorText`, which renders them with the UI
// translation, joined with the UI language's sentence separator (a space in
// English, none in Japanese and Chinese); in English it returns the composed
// English text. Templates
// reach it through `app.config.globalProperties` (app.js).
import { renderErrorText } from './error-normalization.js';
import { itemListText, listText, sentenceSeparator, tm } from './i18n.js';

/** @import { ErrorTextList, ErrorTextPart } from './error-normalization.js' */

/**
 * A list value as the UI language lists shown items.
 * @param {ErrorTextList} list
 * @returns {string}
 */
const shownList = ({ items, and }) => (and ? listText(items) : itemListText(items));

/**
 * The shown text of error text parts.
 * @param {readonly ErrorTextPart[]} parts
 * @returns {string}
 */
export const errorText = (parts) => renderErrorText(parts, tm, sentenceSeparator(), shownList);

/**
 * The shown text of one value of a part: a part through its key, a list as the UI language lists items.
 * @param {ErrorTextPart | ErrorTextList} value
 * @returns {string}
 */
export const errorValueText = (value) => ('items' in value ? shownList(value) : errorText([value]));
