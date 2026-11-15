// @ts-check
// The Web owner of selector text: a record selector (`#N`) and a feature
// selector (`KEY=VALUE`, `hash=<hash>`), as the CLI, the tables and the Session
// read them, and the forms the help and error texts show. Python owns the
// syntax (gbdraw/io/record_select.py, gbdraw/features/selector_values.py) and
// the generated contracts carry it here. Some older writers still write `#`
// or `hash` by hand; they move here when they are next changed. A leaf that
// imports only the contracts, so the error wording owner can read it.
import { FEATURE_HASH_SELECTOR_KEY, RECORD_INDEX_SELECTOR_PREFIX } from '../contracts/schema.generated.js';

/**
 * The record selector of the record at `position` (from 1), or the form of
 * one when `position` is a `<word>`.
 * @param {number | string} position
 * @returns {string}
 */
export const recordIndexSelector = (position) => `${RECORD_INDEX_SELECTOR_PREFIX}${position}`;

/**
 * A `KEY=VALUE` feature selector.
 * @param {string} key
 * @param {string} value
 * @returns {string}
 */
export const featureSelector = (key, value) => `${key}=${value}`;

/**
 * The feature selector that names a feature by its hash.
 * @param {string} hash
 * @returns {string}
 */
export const hashSelector = (hash) => featureSelector(FEATURE_HASH_SELECTOR_KEY, hash);

/**
 * The selector forms that the UI sentences of more than one key show: code
 * text, never translated and never written into a key (OV-2790). A key names
 * the form by a placeholder and its call passes the form as the value, which
 * utils/i18n.js inserts as one run that a right-to-left sentence does not
 * reorder.
 */
export const SELECTOR_TEXT = Object.freeze({
  // The record selector of a record's position among the records.
  recordNumber: recordIndexSelector('<number>'),
  // A feature selector, as an example and as the form that names one feature (Align at feature).
  featureSelectorExample: featureSelector('gene', 'COX1'),
  featureSelectorForm: featureSelector('locus_tag', '<tag>'),
  // The feature selector key of a feature hash.
  hashSelector: hashSelector('')
});
