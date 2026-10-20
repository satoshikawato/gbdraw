// @ts-check
// The one display of user-facing error text (D-17). The wording owner,
// `normalizeUserFacingError` (R6), returns the English `summary` and the parts
// it is composed of (`summaryParts`; `noteParts` of a live-edit failure). A
// display shows the parts through `errorText`, which renders them with the UI
// translation; in English it returns the composed English text. Templates
// reach it through `app.config.globalProperties` (app.js).
import { renderErrorText } from './error-normalization.js';
import { tm } from './i18n.js';

/** @import { ErrorTextPart } from './error-normalization.js' */

/**
 * The shown text of error text parts.
 * @param {readonly ErrorTextPart[]} parts
 * @returns {string}
 */
export const errorText = (parts) => renderErrorText(parts, tm);
