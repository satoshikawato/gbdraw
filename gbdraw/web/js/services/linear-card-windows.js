// @ts-check
// The windows of a Linear record card as its Region block and its Windows list edit them (PD-OI-124): the one owner
// of a card's window list. A card stores one window in its crop fields (`region_*`) and two or more, or one across the
// origin, in `windows` (`createLinearSeq`, `cardCropWindow`); every edit reads the list here, changes it, and stores
// it back in that form, so the card never holds a window in both. An End typed before the Start without "Across the
// origin" is stored in order with Reverse complement ticked (OV-2884): the figure stays the reversed span it was.
import { cardCropWindow } from './drawing/inputs.js';
/** @import { LinearSeq, LinearSeqWindow } from './drawing/inputs.js' */

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
 * Sets fields of the card's window `index` (one History step: the caller's control). False when nothing changed.
 * @param {CardCrop} card
 * @param {number} index
 * @param {CardWindowEdit} edit
 */
export const editCardWindow = (card, index, edit) => {
  const windows = cardWindowList(card);
  const current = windows[index];
  if (!current) return false;
  const next = normalizedCardWindow({ ...current, ...edit });
  if (['start', 'end', 'reverseComplement', 'acrossOrigin'].every(
    (field) => Object.is(next[/** @type {keyof LinearSeqWindow} */ (field)], current[/** @type {keyof LinearSeqWindow} */ (field)])
  )) return false;
  windows[index] = next;
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
 * The bases a window across the origin of a record of `recordLength` draws (Start through the origin to End), or
 * null when it does not cross or the length is unknown.
 * @param {LinearSeqWindow} window
 * @param {number | null} recordLength
 */
export const originWindowLength = (window, recordLength) => (
  window.acrossOrigin && window.start !== null && window.end !== null && window.start > window.end
    && recordLength !== null && window.start <= recordLength
    ? recordLength - window.start + 1 + window.end
    : null
);

/**
 * The bases of the record between two windows of the list that are next to each other and run in the same direction
 * (PD-OI-124): the interval from the earlier window's end to the later one's start in the record's own order, when
 * no other window of the record lies in it. A circular record's interval may run through the origin. `after` is the
 * index of the first window of the pair.
 * @param {readonly LinearSeqWindow[]} windows
 * @param {{ recordLength: number | null, circular: boolean }} record
 * @returns {{ after: number, bp: number }[]}
 */
export const cardWindowGaps = (windows, { recordLength, circular }) => windows.slice(1).flatMap((right, index) => {
  const left = windows[index];
  const spanned = [left, right].every((window) => window.start !== null && window.end !== null);
  if (!spanned || left.reverseComplement !== right.reverseComplement || recordLength === null) return [];
  const [earlier, later] = left.reverseComplement ? [right, left] : [left, right];
  const from = /** @type {number} */ (earlier.end) + 1;
  const to = /** @type {number} */ (later.start) - 1;
  const bp = circular ? (((to - from + 1) % recordLength) + recordLength) % recordLength : to - from + 1;
  if (bp <= 0 || (!circular && to > recordLength)) return [];
  /** @param {number} base */
  const inGap = (base) => (circular ? ((base - from) % recordLength + recordLength) % recordLength < bp : base >= from && base <= to);
  const others = windows.filter((window) => window !== left && window !== right);
  if (others.some((window) => [window.start, window.end].some((base) => base !== null && inGap(base)))) return [];
  return [{ after: index, bp }];
});

/**
 * The card copy Generate draws: each listed window with both ends empty spans the whole record (1 to
 * `recordLength`), as the card says until its Start and End are set (PD-OI-124). A card without such a window is
 * returned as it is; without a length the empty window stays, and the request refuses it.
 * @template {{ windows?: readonly LinearSeqWindow[] }} C
 * @param {C} card
 * @param {number | null} recordLength
 * @returns {C}
 */
export const resolveWholeRecordWindows = (card, recordLength) => {
  const windows = card.windows || [];
  const empty = (/** @type {LinearSeqWindow} */ window) => window.start === null && window.end === null;
  if (recordLength === null || !windows.some(empty)) return card;
  return {
    ...card,
    windows: windows.map((window) => (empty(window) ? { ...window, start: 1, end: recordLength, acrossOrigin: false } : window))
  };
};

/**
 * The length of the record a drawn card reads, from Generate's Linear record catalog (`sourceIndex`, the card's place
 * among the drawn cards), or null when the catalog does not give one.
 * @param {{ records?: readonly { sourceIndex?: unknown, recordLength?: unknown }[] } | null | undefined} catalog
 * @param {number} sourceIndex
 * @returns {number | null}
 */
export const catalogRecordLength = (catalog, sourceIndex) => {
  const record = (catalog?.records || []).find((entry) => Number(entry?.sourceIndex) === sourceIndex);
  const length = Number(record?.recordLength);
  return Number.isSafeInteger(length) && length > 0 ? length : null;
};
