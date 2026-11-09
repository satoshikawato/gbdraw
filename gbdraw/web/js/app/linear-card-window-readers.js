// @ts-check
// The windows of the shown drawing's Linear cards as the card and Record Layout show them (PD-OI-124): each window's
// fields, its circled number in the drawing, the length of a window across the origin, the bases left out between
// two windows, and the Record Layout's cells. The shown-drawing readers pass the state, the shown drawing of a mode,
// and the record a card draws as ports; the edits are services/linear-card-windows.js.
import { cardWindowGaps, cardWindowList, originWindowLength } from '../services/linear-card-windows.js';
import { linearLayoutCells } from '../services/linear-record-layout.js';
import { t } from '../utils/i18n.js';
/** @import { DiagramMode, DrawingState, LinearSeq, WebState } from '../state.js' */
/** @import { LinearSeqWindow } from '../services/drawing/inputs.js' */
/** @import { SequenceRecord } from '../services/record-discovery.js' */

/**
 * One window of a Linear card as its Region block or Windows list shows it (PD-OI-124): the window's fields, its
 * `index` in the card's list, its `number` in the drawing, the length text of a window across the origin, and the
 * text of the bases left out before the next window ('' when none).
 * @typedef {LinearSeqWindow & { index: number, number: number, originText: string, gapText: string }} LinearCardWindowRow
 */
/** @typedef {{ circular: boolean, windows: LinearCardWindowRow[] }} LinearCardWindowView */
/**
 * One cell of the Record Layout panel: the card (`uid`, `cardIndex`), its window (null for a whole card), a stable
 * `key`, and the window's circled number (null for a whole card).
 * @typedef {{ uid: string, cardIndex: number, window: number | null, key: string, number: number | null }} LinearLayoutCellRow
 */

/**
 * @typedef {object} LinearCardWindowReaderPorts
 * @property {WebState} state
 * @property {(drawingMode: DiagramMode) => DrawingState | null} shownDrawingOfMode The shown drawing when it is of that mode.
 * @property {(drawing: DrawingState, seq: LinearSeq) => SequenceRecord | null} linearSeqRecord The record a card draws, once read.
 */

/** @param {LinearCardWindowReaderPorts} ports */
export const createLinearCardWindowReaders = ({ state, shownDrawingOfMode, linearSeqRecord }) => {
  /**
   * The number each window of a card carries in the drawing: its place among the drawn units of its row when the
   * Record Layout arranges the rows, else its place in the card's list (the card is its own row). A card drawn in one
   * window, or not drawn, numbers its list.
   * @param {DrawingState} drawing @param {LinearSeq} seq @param {number} count
   */
  const windowNumbers = (drawing, seq, count) => {
    const cells = count > 1 && drawing.linearRecordLayoutEnabled.value
      ? linearLayoutCells(drawing.drawnLinearSeqs.value, drawing.linearRecordRows) : [];
    return Array.from({ length: count }, (_, index) => (
      cells.find((cell) => cell.uid === seq.uid && cell.window === index)?.position ?? index + 1
    ));
  };
  /**
   * The Region block or the Windows list of a Linear card of the shown drawing (PD-OI-124): each window's fields, its
   * number in the drawing, the length of a window across the origin, and the bases left out between two windows that
   * are next to each other and run in the same direction; `circular`, whether the record is circular (its windows
   * offer "Across the origin").
   * @param {LinearSeq} seq
   * @returns {LinearCardWindowView}
   */
  const linearCardWindowView = (seq) => {
    const drawing = state.activeDrawing();
    const record = linearSeqRecord(drawing, seq);
    const recordLength = typeof record?.recordLength === 'number' ? record.recordLength : null;
    const circular = record?.detectedTopology === 'circular';
    const windows = cardWindowList(seq);
    const numbers = windowNumbers(drawing, seq, windows.length);
    const gaps = cardWindowGaps(windows, { recordLength, circular });
    /** @param {number} value */
    const bp = (value) => value.toLocaleString('en-US');
    return {
      circular,
      windows: windows.map((window, index) => {
        const length = originWindowLength(window, recordLength);
        const gap = gaps.find((entry) => entry.after === index);
        return {
          ...window, index, number: numbers[index],
          originText: length === null || window.start === null || window.end === null || recordLength === null ? ''
            : t('{start} → {recordEnd} → 1 → {end}: {length} bp across the origin.', {
                start: bp(window.start), recordEnd: bp(recordLength), end: bp(window.end), length: bp(length)
              }),
          gapText: gap ? t('{bp} bp of the record lie between windows {first} and {second} and are left out.', {
            bp: bp(gap.bp), first: numbers[index], second: numbers[index + 1]
          }) : ''
        };
      })
    };
  };
  /**
   * The cells of the shown Linear drawing's Record Layout: one per card in list order while no card is drawn in
   * several windows (as before PD-OI-124), else one per window in drawing order. A cell's `number` is the circled
   * number of a window (null for a whole card), counted among the drawn units of its row.
   * @returns {LinearLayoutCellRow[]}
   */
  const linearLayoutCellRows = () => {
    const drawing = shownDrawingOfMode('linear');
    if (!drawing) return [];
    const cards = drawing.linearSeqs;
    const windowed = cards.some((card) => card.windows.length > 1);
    const drawn = windowed ? linearLayoutCells(drawing.drawnLinearSeqs.value, drawing.linearRecordRows) : [];
    const cells = windowed ? linearLayoutCells(cards, drawing.linearRecordRows)
      : cards.map((card, cardIndex) => ({ uid: card.uid, cardIndex, window: null, row: 0, position: 0, placed: false, column: null }));
    return cells.map(({ uid, cardIndex, window }) => ({
      uid, cardIndex, window,
      key: `${uid}:${window ?? ''}`,
      number: window === null ? null
        : drawn.find((cell) => cell.uid === uid && cell.window === window)?.position ?? window + 1
    }));
  };
  return { linearCardWindowView, linearLayoutCellRows };
};
