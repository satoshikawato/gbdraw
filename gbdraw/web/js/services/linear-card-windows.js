// @ts-check
// The windows of a Linear record card as its Region block and its Windows list edit them (PD-OI-124): the one owner
// of a card's window list. A card stores one window in its crop fields (`region_*`) and two or more, or one across the
// origin, in `windows` (`createLinearSeq`, `cardCropWindow`); every edit reads the list here, changes it, and stores
// it back in that form, so the card never holds a window in both. An End typed before the Start without "Across the
// origin" is stored in order with Reverse complement ticked (OV-2884): the figure stays the reversed span it was.
import { cardCropWindow } from './drawing/inputs.js';
import { cropFormReverseComplement, linearCardCropForm, orientedCropForm } from './record-display-model.js';
import { cloneJsonData } from './json-clone.js';
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
 * A window as the card stores it: without "Across the origin", an End before the Start is the reversed span, stored
 * in order with Reverse complement ticked (OV-2884, D-55); with it, the span is kept as typed.
 * @param {LinearSeqWindow} window
 * @returns {LinearSeqWindow}
 */
export const normalizedCardWindow = (window) => (
  !window.acrossOrigin && window.start !== null && window.end !== null && window.start > window.end
    ? { ...window, start: window.end, end: window.start, reverseComplement: true }
    : window
);

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
 * What the page knows of the record a card draws: its length (null until read) and its ID ('' until read).
 * @typedef {{ length: number | null, id: string }} CardRecordFacts
 */
/**
 * The records of Generate's Linear record catalog, as this module reads them.
 * @typedef {{ records?: readonly { sourceIndex?: unknown, recordLength?: unknown, recordId?: unknown }[] } | null | undefined} LinearRecordCatalogRecords
 */
/** A record not read yet. @type {CardRecordFacts} */
export const UNKNOWN_CARD_RECORD = Object.freeze({ length: null, id: '' });

/**
 * The record a drawn card reads, from Generate's Linear record catalog (`sourceIndex`, the card's place among the
 * drawn cards).
 * @param {LinearRecordCatalogRecords} catalog
 * @param {number} sourceIndex
 * @returns {CardRecordFacts}
 */
export const catalogCardRecord = (catalog, sourceIndex) => {
  const record = (catalog?.records || []).find((entry) => Number(entry?.sourceIndex) === sourceIndex);
  return sequenceCardRecord(record);
};

/**
 * The facts of a record as a record list or the catalog holds it.
 * @param {{ recordLength?: unknown, recordId?: unknown } | null | undefined} record
 * @returns {CardRecordFacts}
 */
export const sequenceCardRecord = (record) => {
  const length = Number(record?.recordLength);
  return {
    length: Number.isSafeInteger(length) && length > 0 ? length : null,
    id: typeof record?.recordId === 'string' ? record.recordId : ''
  };
};
