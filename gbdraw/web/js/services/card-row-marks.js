// @ts-check
// The mark a row of each file card carries, and the removal of it. One owner:
// the card helpers (services/file-imports.js) and the label-tables slice
// (an edit of a Whitelist row) both read it. It imports nothing, so a slice may.
/** @import { CardMarkedRow, RowFileCard } from './file-imports.js' */

/**
 * The mark a row of each card's file carries.
 * @typedef {{ t_color: { fromFile: true }, whitelist: { importedFrom: 'whitelist' }, qualifier_priority: { importedFrom: 'qualifier_priority' } }} CardRowMarks
 */
/** @type {Readonly<CardRowMarks>} */
export const CARD_ROW_MARKS = Object.freeze({
  t_color: Object.freeze({ fromFile: true }),
  whitelist: Object.freeze({ importedFrom: 'whitelist' }),
  qualifier_priority: Object.freeze({ importedFrom: 'qualifier_priority' })
});

/**
 * Makes `row` the user's: it loses `card`'s mark (an edit of the row).
 * @param {CardMarkedRow} row Changed in place.
 * @param {RowFileCard} card
 */
export const unmarkCardRow = (row, card) => {
  Object.keys(CARD_ROW_MARKS[card]).forEach((key) => { delete row[/** @type {keyof CardMarkedRow} */ (key)]; });
};
