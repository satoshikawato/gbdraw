// @ts-check
// A list cut to a line limit, and the one line that stands for the rest (OV-3495): the Load repair dialog and the
// preview's input-summary notice show long lists this way.
import { tn } from './i18n.js';

/**
 * @template T
 * @typedef {{ readonly shown: readonly T[], readonly rest: readonly T[] }} LimitedList
 *   A list cut to its line limit: the items it shows, and the rest that one `moreLine` stands for.
 */

/**
 * A list within `limit` lines: all of it when it fits, else its first `limit - 1` items and the rest, which one
 * `moreLine(rest.length)` stands for. `limit` is a whole number of at least 1 (1 shows only the `moreLine`).
 * @template T
 * @param {readonly T[]} items
 * @param {number} limit
 * @returns {LimitedList<T>}
 */
export const limitedList = (items, limit) => Object.freeze(items.length > limit
  ? { shown: Object.freeze(items.slice(0, limit - 1)), rest: Object.freeze(items.slice(limit - 1)) }
  : { shown: Object.freeze([...items]), rest: Object.freeze([]) });

/**
 * The line that stands for the rest of a `limitedList`, in the shown UI language.
 * @param {number} count `rest.length`
 * @returns {string}
 */
export const moreLine = (count) => tn(count, 'and one more', 'and {count} more');
