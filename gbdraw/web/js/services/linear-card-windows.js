// @ts-check
// The windows of a Linear record card as its Region block and its Windows list edit them (PD-OI-124): the one owner
// of a card's window list. A card stores one window in its crop fields (`region_*`) and two or more, or one across the
// origin, in `windows` (`createLinearSeq`, `cardCropWindow`); every edit reads the list here, changes it, and stores
// it back in that form, so the card never holds a window in both. An End typed before the Start without "Across the
// origin" is stored in order with Reverse complement ticked (OV-2884): the figure stays the reversed span it was.
import { cardCropWindow } from './drawing/inputs.js';
import { cropFormReverseComplement, linearCardCropForm, orientedCropForm } from './record-display-model.js';
import { cloneJsonData } from './json-clone.js';
import { resolveDisambiguatedRecordSelection } from './record-options.js';
/** @import { LinearSeq, LinearSeqWindow } from './drawing/inputs.js' */
/** @import { LinearCardOrientation } from './record-display-model.js' */

/** The crop fields of a card that this module reads and writes. @typedef {Pick<LinearSeq, 'region_start' | 'region_end' | 'region_reverse' | 'windows'>} CardCrop */
/** The fields of a window an edit may set. @typedef {Partial<Pick<LinearSeqWindow, 'start' | 'end' | 'reverseComplement' | 'acrossOrigin'>>} CardWindowEdit */

/** @param {unknown} value @returns {number | null} */
const endpoint = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  return Number.isFinite(numeric) ? numeric : null;
};

/**
 * The card's windows in drawing order: its list, or its crop fields as one window (both ends empty when the card
 * draws the whole record).
 * @param {CardCrop} card
 * @returns {LinearSeqWindow[]}
 */
export const cardWindowList = (card) => (card.windows.length
  ? card.windows.map((window) => ({ ...window }))
  : [{
      start: endpoint(card.region_start), end: endpoint(card.region_end),
      reverseComplement: Boolean(card.region_reverse), column: null, acrossOrigin: false
    }]);

/**
 * Stores `windows` on the card: one window that does not cross the origin in its crop fields, else the list.
 * @param {CardCrop} card
 * @param {LinearSeqWindow[]} windows
 */
const storeCardWindows = (card, windows) => {
  const only = cardCropWindow(windows);
  card.region_start = only ? only.start : null;
  card.region_end = only ? only.end : null;
  card.region_reverse = only ? only.reverseComplement : false;
  card.windows = only ? [] : windows;
};

/**
 * The span a window's Start, End and "Across the origin" give (PD-OI-124): with the box, a Start greater than the End
 * runs from the Start through the origin to the End, and a Start not greater than the End is no window (null: Generate
 * and Find in other files refuse it); without it, an End before the Start reverses the span. The one owner of that
 * rule for the card, the request and the Find in other files reference rows.
 * @param {number} start @param {number} end @param {boolean} acrossOrigin
 * @returns {{ start: number, end: number, reverseComplement: boolean } | null}
 */
export const windowSpanOf = (start, end, acrossOrigin) => (acrossOrigin
  ? (start > end ? { start, end, reverseComplement: false } : null)
  : { start: Math.min(start, end), end: Math.max(start, end), reverseComplement: start > end });

/**
 * A window as the card stores it: without "Across the origin", an End before the Start is the reversed span, stored
 * in order with Reverse complement ticked (OV-2884, D-55, `windowSpanOf`); with it, the span is kept as typed.
 * @param {LinearSeqWindow} window
 * @returns {LinearSeqWindow}
 */
export const normalizedCardWindow = (window) => {
  if (window.acrossOrigin || window.start === null || window.end === null) return window;
  const span = /** @type {{ start: number, end: number, reverseComplement: boolean }} */ (windowSpanOf(window.start, window.end, false));
  return span.reverseComplement ? { ...window, start: span.start, end: span.end, reverseComplement: true } : window;
};

/**
 * What an edit of the card's window `index` changes: 'crop' when the window's span changes (Start, End, Across the
 * origin, or an End typed before the Start, which also ticks Reverse complement), 'orientation' when only its Reverse
 * complement changes (the direction an alignment records, as Align and Reset write it on the card), null when nothing
 * changes.
 * @param {CardCrop} card
 * @param {number} index
 * @param {CardWindowEdit} edit
 * @returns {'crop' | 'orientation' | null}
 */
export const cardWindowEditKind = (card, index, edit) => {
  const current = cardWindowList(card)[index];
  if (!current) return null;
  const next = normalizedCardWindow({ ...current, ...edit });
  /** @param {'start' | 'end' | 'reverseComplement' | 'acrossOrigin'} field */
  const changed = (field) => !Object.is(next[field], current[field]);
  if (changed('start') || changed('end') || changed('acrossOrigin')) return 'crop';
  return changed('reverseComplement') ? 'orientation' : null;
};

/**
 * Sets fields of the card's window `index` (one History step: the caller's control). False when nothing changed.
 * @param {CardCrop} card
 * @param {number} index
 * @param {CardWindowEdit} edit
 */
export const editCardWindow = (card, index, edit) => {
  if (!cardWindowEditKind(card, index, edit)) return false;
  const windows = cardWindowList(card);
  windows[index] = normalizedCardWindow({ ...windows[index], ...edit });
  storeCardWindows(card, windows);
  return true;
};

/**
 * Adds a window after the card's last one, with Start and End empty (the whole record until they are set) and
 * `column`, its place in the arranged row (null when the rows are not arranged).
 * @param {CardCrop} card
 * @param {number | null} column
 */
export const addCardWindow = (card, column) => {
  storeCardWindows(card, [...cardWindowList(card),
    { start: null, end: null, reverseComplement: false, column, acrossOrigin: false }]);
};

/**
 * Removes the card's window `index`; the last window stays. False when nothing was removed.
 * @param {CardCrop} card
 * @param {number} index
 */
export const removeCardWindow = (card, index) => {
  const windows = cardWindowList(card);
  if (windows.length < 2 || !windows[index]) return false;
  windows.splice(index, 1);
  storeCardWindows(card, windows);
  return true;
};

/**
 * Puts back the orientation fields `captureLinearCardOrientation` (services/record-display-model.js) kept, as they
 * were.
 * @template {{ reverseComplement: boolean }} W
 * @param {LinearCardOrientation<W>} card
 * @param {LinearCardOrientation<W>} saved
 */
export const restoreLinearCardOrientation = (card, saved) => {
  card.windows = cloneJsonData(saved.windows);
  card.region_reverse = saved.region_reverse;
};

/**
 * Draws a Linear File card's record in the absolute direction `reverseComplement` where the card holds its
 * orientation (`orientedCropForm`): on the windows it lists, else on its Reverse complement (Record options Reverse,
 * Align's flip, the record popup's Reverse complement). A card drawn so already is left as it is; false then.
 * @template {{ reverseComplement: boolean }} W
 * @param {LinearCardOrientation<W>} card
 * @param {boolean} reverseComplement
 */
export const orientLinearCard = (card, reverseComplement) => {
  const form = linearCardCropForm(card);
  if (cropFormReverseComplement(form) === reverseComplement) return false;
  const oriented = orientedCropForm(form, reverseComplement);
  card.windows = oriented.windows;
  card.region_reverse = oriented.reverse;
  return true;
};

/**
 * A listed window as Generate draws it: with Start and End both empty it spans the whole record, 1 to `recordLength`,
 * as the card says until they are set (PD-OI-124); any other window, or an empty one whose record length is not known
 * yet, is returned as it is (the request then refuses it, naming the card). The one owner of that fill: every reader
 * of a card's windows reads them through `linearWindowsPayload` (services/session-request.js), which calls it.
 * @template {Pick<LinearSeqWindow, 'start' | 'end'> & { acrossOrigin?: boolean }} W
 * @param {W} window
 * @param {number | null} recordLength
 * @returns {W}
 */
export const wholeRecordWindow = (window, recordLength) => (
  window.start === null && window.end === null && recordLength !== null
    ? { ...window, start: 1, end: recordLength, acrossOrigin: false }
    : window
);

/**
 * What the page knows of the records a card draws: `count`, how many (one for a chosen record or a file of one
 * record, every record of the file otherwise; null until read), and the `length` and `id` of the one record it draws
 * (null and '' until read, or while it draws several). `fileCount` is how many records the card's file holds (null
 * until read): the request names the record of a file of one record as the CLI does, by no selector
 * (`cliRecordSelector` in session-request.js).
 * @typedef {{ length: number | null, id: string, count: number | null, fileCount: number | null }} CardRecordFacts
 */
/**
 * The records of Generate's Linear record catalog, as this module reads them: each with the number of records its
 * card's file holds (`sourceRecordCount`).
 * @typedef {{ records?: readonly { sourceIndex?: unknown, recordLength?: unknown, recordId?: unknown,
 *   sourceRecordCount?: unknown }[] } | null | undefined} LinearRecordCatalogRecords
 */
/** A record not read yet. @type {CardRecordFacts} */
export const UNKNOWN_CARD_RECORD = Object.freeze({ length: null, id: '', count: null, fileCount: null });

/** @param {unknown} value @returns {number | null} */
const positiveCount = (value) => {
  const count = Number(value);
  return Number.isSafeInteger(count) && count > 0 ? count : null;
};

/**
 * The facts of the records a card draws, as a record list or the catalog holds them: the one owner for Generate
 * (`catalogCardRecord`) and the card while editing (the record selector's records, `cardFileRecordFacts`).
 * @param {readonly ({ recordLength?: unknown, recordId?: unknown } | null | undefined)[]} records
 * @param {number | null} fileCount The records of the card's file, null until read.
 * @returns {CardRecordFacts}
 */
export const cardRecordFacts = (records, fileCount) => {
  if (records.length !== 1) return { ...UNKNOWN_CARD_RECORD, count: records.length || null, fileCount };
  return {
    length: positiveCount(records[0]?.recordLength),
    id: typeof records[0]?.recordId === 'string' ? records[0].recordId : '',
    count: 1,
    fileCount
  };
};

/**
 * The records a card draws of its file's `records`: the one its record selector names, else every record when it
 * names none ([] when it names no single record of the file).
 * @template {{ recordId?: unknown, selector?: unknown }} R
 * @param {readonly R[]} records
 * @param {unknown} selector The card's `region_record_id`.
 * @returns {readonly R[]}
 */
export const cardDrawnRecords = (records, selector) => {
  const selection = resolveDisambiguatedRecordSelection(records, selector);
  return selection.record ? [selection.record] : selection.status === 'unspecified' ? records : [];
};

/**
 * The facts of a card from its file's records as read (`[]` until read): the card while editing and a request
 * rebuilt outside the app (a Gallery publication) read them so.
 * @param {readonly { recordId?: unknown, selector?: unknown, recordLength?: unknown }[]} records
 * @param {unknown} selector The card's `region_record_id`.
 * @returns {CardRecordFacts}
 */
export const cardFileRecordFacts = (records, selector) => cardRecordFacts(cardDrawnRecords(records, selector),
  records.length || null);

/**
 * The records a drawn card reads, from Generate's Linear record catalog (`sourceIndex`, the card's place among the
 * drawn cards; the catalog holds the records its selector chose).
 * @param {LinearRecordCatalogRecords} catalog
 * @param {number} sourceIndex
 * @returns {CardRecordFacts}
 */
export const catalogCardRecord = (catalog, sourceIndex) => {
  const records = (catalog?.records || []).filter((entry) => Number(entry?.sourceIndex) === sourceIndex);
  return cardRecordFacts(records, positiveCount(records[0]?.sourceRecordCount));
};
