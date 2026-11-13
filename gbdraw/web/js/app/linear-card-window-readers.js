// @ts-check
// The windows of the shown drawing's Linear cards as the card and Record Layout show them (PD-OI-124): each window's
// fields, its circled number in the drawing, the length of a window across the origin, the bases left out between
// two windows, and the Record Layout's cells. The numbers and the omitted lengths are Python's (the Worker helper
// `planCardWindows` reads the records as the request writes them, with the drawing's row order and captions); until
// it answers, the card shows neither. The shown-drawing readers pass the state, the shown drawing of a mode, the
// record a card draws, the helper and the record display rows as ports; the edits are services/linear-card-windows.js.
import { cardWindowList } from '../services/linear-card-windows.js';
import { windowBases } from '../services/record-options.js';
import { linearLayoutCells } from '../services/linear-record-layout.js';
import { numberText, t } from '../utils/i18n.js';
/** @import { DiagramMode, DrawingState, LinearSeq, WebState } from '../state.js' */
/** @import { LinearSeqWindow } from '../services/drawing/inputs.js' */
/** @import { CardRecordFacts } from '../services/linear-card-windows.js' */
/** @import { SequenceRecord } from '../services/record-discovery.js' */
/** @import { RequestRecordDisplayRow, linearCardWindowPlanRecords } from '../services/session-request.js' */
/** @import { CardWindowPlan, CardWindowRecordInput } from '../contracts/card-windows.generated.js' */

const { computed, ref, watch } = window.Vue;

/**
 * One window of a Linear card as its Region block or Windows list shows it (PD-OI-124): the window's fields, its
 * `index` in the card's list, its circled `number` in the drawing (null until Python answers, or when the drawing
 * does not number it), `label`, the number its controls are named by (its circled number, else its place in the
 * list), the length text of a window across the origin, and the text of the bases left out before the window drawn
 * after it ('' when none).
 * @typedef {LinearSeqWindow & { index: number, number: number | null, label: number, originText: string, gapText: string }} LinearCardWindowRow
 */
/** @typedef {{ circular: boolean, windows: LinearCardWindowRow[] }} LinearCardWindowView */
/**
 * One cell of the Record Layout panel: the card (`uid`, `cardIndex`), its window (null for a whole card), a stable
 * `key`, the window's circled number (null for a whole card, or until Python answers) and `label`, the number its
 * Row field is named by.
 * @typedef {{ uid: string, cardIndex: number, window: number | null, key: string, number: number | null, label: number | null }} LinearLayoutCellRow
 */

/**
 * @typedef {object} LinearCardWindowReaderPorts
 * @property {WebState} state
 * @property {(drawingMode: DiagramMode) => DrawingState | null} shownDrawingOfMode The shown drawing when it is of that mode.
 * @property {(drawing: DrawingState, seq: LinearSeq) => SequenceRecord | null} linearSeqRecord The record a card draws, once read.
 * @property {(drawing: DrawingState, seq: LinearSeq) => CardRecordFacts} linearCardRecord What a card's record list says
 *   of the records it draws (app/shown-drawing-readers.js), the one owner of a card's record ID and length.
 * @property {(payload: { records: CardWindowRecordInput[] }) => Promise<{ result?: unknown }>} planCardWindows
 *   The Worker helper that plans the windows' numbers and omitted lengths.
 * @property {() => RequestRecordDisplayRow[]} recordDisplayRows The record display rows the request reads.
 * @property {typeof linearCardWindowPlanRecords} cardWindowPlanRecords The records the helper reads, as the request
 *   writes them (`linearCardWindowPlanRecords` of services/session-request.js).
 */

/** @param {unknown} result @returns {CardWindowPlan | null} */
const cardWindowPlanOf = (result) => {
  const plan = /** @type {Partial<CardWindowPlan> | null} */ (result && typeof result === 'object' ? result : null);
  return plan && Array.isArray(plan.windows) && Array.isArray(plan.gaps) ? /** @type {CardWindowPlan} */ (plan) : null;
};

/** @param {LinearCardWindowReaderPorts} ports */
export const createLinearCardWindowReaders = ({
  state, shownDrawingOfMode, linearSeqRecord, linearCardRecord, planCardWindows, cardWindowPlanRecords,
  recordDisplayRows
}) => {
  /** @param {DrawingState} drawing @param {LinearSeq} seq */
  const cardRecordOf = (drawing, seq) => ({
    ...linearCardRecord(drawing, seq),
    circular: linearSeqRecord(drawing, seq)?.detectedTopology === 'circular'
  });
  // The payload Python plans for the shown Linear drawing ('' while no drawn card has several windows, a record is not
  // read, or a window is one Generate refuses).
  const planKey = computed(() => {
    const drawing = shownDrawingOfMode('linear');
    const seqs = drawing ? drawing.drawnLinearSeqs.value : [];
    if (!drawing || !seqs.some((seq) => seq.windows.length > 1)) return '';
    const records = cardWindowPlanRecords(drawing, seqs, (seq) => cardRecordOf(drawing, seq), recordDisplayRows());
    return records ? JSON.stringify(records) : '';
  });
  /** @type {{ value: { key: string, plan: CardWindowPlan | null } | null }} */
  const answered = ref(null);
  watch(() => planKey.value, /** @param {string} key */ async (key) => {
    if (!key) return;
    let plan = null;
    try {
      plan = cardWindowPlanOf((await planCardWindows({ records: JSON.parse(key) })).result);
    } catch {
      plan = null;
    }
    if (planKey.value === key) answered.value = { key, plan };
  }, { immediate: true });
  /** Python's plan of the shown payload, or null until it answers. */
  const shownPlan = () => {
    const key = planKey.value;
    const reply = answered.value;
    return key && reply && reply.key === key ? reply.plan : null;
  };
  /** @param {CardWindowPlan | null} plan @param {string} uid @param {number} window */
  const numberOf = (plan, uid, window) => plan?.windows.find((entry) => entry.recordKey === uid && entry.window === window)?.number ?? null;

  /**
   * The Region block or the Windows list of a Linear card of the shown drawing (PD-OI-124): each window's fields, its
   * number in the drawing, the length of a window across the origin, and the bases left out between two windows that
   * stand next to each other in the drawing and run in the same direction; `circular`, whether the record is circular
   * (its windows offer "Across the origin").
   * @param {LinearSeq} seq
   * @returns {LinearCardWindowView}
   */
  const linearCardWindowView = (seq) => {
    const drawing = state.activeDrawing();
    const { length: recordLength, circular } = cardRecordOf(drawing, seq);
    const plan = shownPlan();
    const windows = cardWindowList(seq);
    return {
      circular,
      windows: windows.map((window, index) => {
        // A window across the origin: its length runs from its Start through the origin to its End.
        const length = window.acrossOrigin && window.start !== null && window.end !== null && window.start > window.end
          && recordLength !== null && window.start <= recordLength ? windowBases(window.start, window.end, recordLength) : null;
        const number = numberOf(plan, seq.uid, index);
        const gap = plan?.gaps.find((entry) => entry.recordKey === seq.uid && entry.leftWindow === index);
        const second = gap ? numberOf(plan, seq.uid, gap.rightWindow) : null;
        return {
          ...window, index, number, label: number ?? index + 1,
          originText: length === null || window.start === null || window.end === null || recordLength === null ? ''
            : t('{start} → {recordEnd} → 1 → {end}: {length} bp across the origin.', {
                start: numberText(window.start), recordEnd: numberText(recordLength), end: numberText(window.end), length: numberText(length)
              }),
          gapText: gap && number !== null && second !== null
            ? t('{bp} bp of the record lie between windows {first} and {second} and are left out.', {
                bp: numberText(gap.omittedBp), first: number, second
              }) : ''
        };
      })
    };
  };
  /**
   * The cells of the shown Linear drawing's Record Layout: one per card in list order while no card is drawn in
   * several windows (as before PD-OI-124), else one per window in the order the layout places them. A cell's `number`
   * is the circled number of a window in the drawing (null for a whole card, or until Python answers).
   * @returns {LinearLayoutCellRow[]}
   */
  const linearLayoutCellRows = () => {
    const drawing = shownDrawingOfMode('linear');
    if (!drawing) return [];
    const cards = drawing.linearSeqs;
    const windowed = cards.some((card) => card.windows.length > 1);
    const plan = windowed ? shownPlan() : null;
    const cells = windowed ? linearLayoutCells(cards, drawing.linearRecordRows)
      : cards.map((card, cardIndex) => ({ uid: card.uid, cardIndex, window: null }));
    return cells.map(({ uid, cardIndex, window }) => {
      const number = window === null ? null : numberOf(plan, uid, window);
      return {
        uid, cardIndex, window, key: `${uid}:${window ?? ''}`, number,
        label: window === null ? null : number ?? window + 1
      };
    });
  };
  return { linearCardWindowView, linearLayoutCellRows };
};
