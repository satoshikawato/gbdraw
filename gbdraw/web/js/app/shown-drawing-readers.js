// @ts-check
// Readers of the shown drawing and of the Result it shows: the record a Linear
// card draws and the notes beside the Region fields, the shown Result's input
// summary, the records the next Linear Generate draws and those a loaded
// Result drew, the wait for a deferred Linear record read, the Web-only moves
// of the canvas, the committed request of the Result a drawing shows, the
// record IDs that Result names its records by, and Wrap Rows as that Result
// drew it. The composition root passes the state, the
// shown drawing of a mode, the Linear record selector's members, and the
// committed request and its Linear record cards as ports.
import { inputSummaryOf } from '../services/input-summary.js';
import { readRegionAxisRecords } from '../services/region-axis.js';
import { cardHoldsRecordKey } from '../services/record-draw-selection.js';
import { requestLinearWrapLeaf, shownLinearWrap, wrapFormValue } from '../services/linear-wrap-form.js';
import {
  circularDrawingRegionNote, circularDrawingWindowRows, linearCardRegion, regionClampNote
} from '../services/record-options.js';
import { createLinearCardWindowReaders } from './linear-card-window-readers.js';
import { cardDrawnRecords, cardFileRecordFacts } from '../services/linear-card-windows.js';
/** @import { DiagramMode, DrawingState, LinearSeq, SessionBusy, WebState } from '../state.js' */
/** @import { CanonicalRenderRequest } from '../contracts/request.generated.js' */
/** @import { LinearDraftRecord } from '../services/linear-wrap-form.js' */
/** @import { SequenceRecord } from '../services/record-discovery.js' */
/** @import { RecordSelectorStatus, RefreshOutcome } from './linear-record-selector.js' */
/** @import { LinearCardWindowReaderPorts } from './linear-card-window-readers.js' */
/** @import { linearDrawingUnitPlaces, linearDrawingWindowsOf } from '../services/session-request.js' */

const { computed } = window.Vue;

/**
 * @typedef {object} ShownDrawingReaderPorts
 * @property {WebState} state
 * @property {(drawingMode: DiagramMode) => DrawingState | null} shownDrawingOfMode The shown drawing when it is of that mode.
 * @property {(drawing: DrawingState, seq: LinearSeq) => SequenceRecord[]} linearRecordsFor
 *   The records the record selector read for a Linear card.
 * @property {(drawing: DrawingState, seq: LinearSeq) => RecordSelectorStatus} linearRecordStatusFor
 * @property {(drawing: DrawingState) => SessionBusy | Promise<RefreshOutcome>} refreshLinearRecords
 *   Reads the records of the drawing's Linear cards.
 * @property {(drawing: DrawingState) => { record: CanonicalRenderRequest['records'][number], seq: LinearSeq }[] | null} committedLinearRecordCards
 *   Each record of the committed Linear request with the card of the drawing that still draws it.
 * @property {() => CanonicalRenderRequest | null} getCommittedCanonicalRenderRequest
 * @property {typeof linearDrawingWindowsOf} windowsOf The windows a Linear card writes in its request
 *   (`linearDrawingWindowsOf` of services/session-request.js).
 * @property {typeof linearDrawingUnitPlaces} unitPlacesOf The places a Linear drawing's request writes for its drawn
 *   cards, which its drawn-unit plan reads (`linearDrawingUnitPlaces` of services/session-request.js).
 * @property {LinearCardWindowReaderPorts['cardWindowPlanRecords']} cardWindowPlanRecords
 * @property {LinearCardWindowReaderPorts['planCardWindows']} planCardWindows
 *   The Worker helper that plans the numbers and omitted lengths of a Linear card's windows.
 * @property {LinearCardWindowReaderPorts['recordDisplayRows']} recordDisplayRows The record display rows the request reads.
 */

/**
 * The record ID a shown Result names a record by, and whether that ID names the
 * record without its `#k` as far as the committed request tells (`idSuffices`):
 * the request selects the record by that ID (the Web writes the record's `#k`
 * instead when the ID repeats in its file or is unsafe), or draws the record's
 * whole file (a card the CLI or the Python API wrote: `all`, with no selector
 * and no window), so the Result names every record of the file and an ID that
 * repeats in it repeats among them.
 * @typedef {{ recordId: string, idSuffices: boolean }} ShownRecordName
 */

/**
 * The record ID a Result names each of its records by, by record key (OV-2959):
 * the one Python drew the record with, as the Result's run metadata keeps it
 * (its track-slot geometry, which a Session keeps), else, for a Result saved
 * before Session 51 (whose run metadata keys no record), the record ID its
 * committed request selects the record by.
 * @param {CanonicalRenderRequest | null} request The committed request of the Result.
 * @param {unknown} trackSlotGeometry The run metadata `trackSlotGeometry` of the Result.
 * @param {number} resultCount
 * @returns {ReadonlyMap<string, ShownRecordName>}
 */
const recordIdsOfResult = (request, trackSlotGeometry, resultCount) => {
  /** @type {Map<string, string>} */
  const selected = new Map();
  /** @type {string[]} */
  const wholeFiles = [];
  for (const { recordKey, cardinality, selector, windows } of request?.records ?? []) {
    if (selector?.kind === 'recordId' && selector.value.trim()) selected.set(recordKey, selector.value.trim());
    if (cardinality === 'all' && !selector && windows.length === 0) wholeFiles.push(recordKey);
  }
  /** @type {Map<string, ShownRecordName>} */
  const ids = new Map([...selected].map(([recordKey, recordId]) => [recordKey, { recordId, idSuffices: true }]));
  for (let index = 0; index < resultCount; index += 1) {
    for (const { recordKey, recordId } of readRegionAxisRecords(trackSlotGeometry, index)) {
      if (recordId) {
        const idSuffices = selected.get(recordKey) === recordId || wholeFiles.some((card) => cardHoldsRecordKey(card, recordKey));
        ids.set(recordKey, { recordId, idSuffices });
      }
    }
  }
  return ids;
};

/**
 * Registers the computeds (`circularRegionNote`, `circularWindowRows`, `shownInputSummary`,
 * `shownLinearWrapReading`); the root calls it where they stood.
 * @param {ShownDrawingReaderPorts} ports
 */
export const createShownDrawingReaders = ({
  state, shownDrawingOfMode, linearRecordsFor, linearRecordStatusFor, refreshLinearRecords, committedLinearRecordCards,
  getCommittedCanonicalRenderRequest, windowsOf, unitPlacesOf, planCardWindows, cardWindowPlanRecords, recordDisplayRows
}) => {
  const { inputSummaries, selectedResultIndex, results } = state;
  /**
   * The records a Linear entry draws, once read: the one its selector names, else every record of its file ([]
   * until read, or when its selector names none).
   * @param {DrawingState} drawing the drawing that holds the card
   * @param {LinearSeq} seq
   */
  const linearSeqDrawnRecords = (drawing, seq) => cardDrawnRecords(linearRecordsFor(drawing, seq), seq.region_record_id);
  /** The one record a Linear entry draws, once read. @param {DrawingState} drawing @param {LinearSeq} seq */
  const linearSeqRecord = (drawing, seq) => (linearSeqDrawnRecords(drawing, seq).length === 1
    ? linearSeqDrawnRecords(drawing, seq)[0] : null);
  /**
   * What a Linear card's record list says of the records it draws (`cardRecordFacts`): the one owner of a card's
   * record ID and length while editing and in Generate, which reads the list after it awaits the read (OV-3227).
   * @param {DrawingState} drawing @param {LinearSeq} seq
   */
  const linearCardRecord = (drawing, seq) => cardFileRecordFacts(linearRecordsFor(drawing, seq), seq.region_record_id);
  // The note beside a Linear entry's Region fields (`regionClampNote`). The
  // request draws the larger coordinate as the region's end (reverse when typed first).
  /** @param {LinearSeq} seq a card of the shown drawing */
  const linearRegionNote = (seq) => {
    const region = linearCardRegion(seq, 1, linearSeqRecord(state.activeDrawing(), seq)?.recordLength);
    return region && !('error' in region) ? regionClampNote(region.drawn) : '';
  };
  // The windows a card of the drawing writes in its request, its record as the record selector read it
  // (`linearDrawingWindowsOf`): the wrap estimate reads them here.
  /** @param {DrawingState} drawing */
  const linearCardWindowsOf = (drawing) => windowsOf(drawing, (card) => linearCardRecord(drawing, card));
  // The places the drawn cards of the drawing write in its request (`linearDrawingUnitPlaces`): the LOSAT job
  // estimate's drawn-unit plan reads them here.
  /** @param {DrawingState} drawing */
  const linearUnitPlacesOf = (drawing) => unitPlacesOf(drawing, (card) => linearCardRecord(drawing, card));
  const { linearCardWindowView, linearLayoutCellRows } = createLinearCardWindowReaders({
    state, shownDrawingOfMode, linearSeqRecord, linearCardRecord, planCardWindows, cardWindowPlanRecords, recordDisplayRows
  });
  const circularRegionNote = computed(() => {
    const drawing = shownDrawingOfMode('circular');
    return drawing ? circularDrawingRegionNote(drawing.circularRecordList.value, drawing) : '';
  });
  // The read-only list of the window across the origin the shown Circular drawing's record keeps (D-45).
  const circularWindowRows = computed(() => {
    const drawing = shownDrawingOfMode('circular');
    return drawing ? circularDrawingWindowRows(drawing.circularRecordList.value, drawing) : [];
  });
  // The input summary of the shown Result: its warnings in the status area, its lines in Run info, its
  // wrapped rows.
  const shownInputSummary = computed(() => inputSummaryOf(
    inputSummaries.value, selectedResultIndex.value, results.value[selectedResultIndex.value]?.name
  ));
  // The records the next Linear Generate draws, as their cards ask for them:
  // each source record's length and the window its card writes in the request
  // (a window across the origin starts after it ends, as Python's summary names
  // it). Known once every drawn record is read and every region can be drawn,
  // else none; a card with several windows cannot be wrapped, so it is unknown.
  /** @param {DrawingState | null | undefined} drawing @returns {LinearDraftRecord[]} */
  const linearDraftRecordsOf = (drawing) => {
    if (!drawing || drawing.mode !== 'linear') return [];
    const cardWindows = linearCardWindowsOf(drawing);
    const records = drawing.drawnLinearSeqs.value.map((seq, index) => {
      const sourceLength = Number(linearSeqRecord(drawing, seq)?.recordLength);
      const region = linearCardRegion(seq, index + 1);
      if (!(sourceLength > 0) || (region !== null && 'error' in region)) return null;
      /** @type {CanonicalRenderRequest['records'][number]['windows']} */
      let windows;
      try {
        windows = cardWindows(seq);
      } catch {
        return null; // A window Generate refuses (REGION_INVALID): nothing is known.
      }
      const [only = null] = windows;
      return windows.length < 2 ? { sourceLength, region: only && { start: only.start, end: only.end } } : null;
    });
    return records.length && records.every((record) => record !== null) ? /** @type {LinearDraftRecord[]} */ (records) : [];
  };
  // Their lengths, while none is cropped: the Web predicts no cropped length.
  /** @param {DrawingState | null | undefined} drawing @returns {number[]} */
  const linearRecordLengthsOf = (drawing) => {
    const records = linearDraftRecordsOf(drawing);
    return records.every((record) => record.region === null) ? records.map((record) => record.sourceLength) : [];
  };
  const drawnLinearRecordLengths = () => linearRecordLengthsOf(state.activeDrawing());
  // The records a shown Result without an input summary (a loaded one) drew:
  // each committed record's region as its request asked for it, with the
  // length of the source record its card reads, matched by record key, so a
  // card edited since Load changes neither (D-73). [] while a record has no
  // card that still draws it, or its length is unknown.
  /** @param {DrawingState} drawing @returns {LinearDraftRecord[]} */
  const committedLinearRecordsOf = (drawing) => {
    const records = (committedLinearRecordCards(drawing) ?? []).map(({ record, seq }) => {
      const sourceLength = Number(linearSeqRecord(drawing, seq)?.recordLength);
      // A record drawn in several windows cannot be wrapped (Python refuses it), so it has no wrap estimate.
      const [region = null] = record.windows.length === 1 ? record.windows : [];
      return sourceLength > 0 && record.windows.length < 2
        ? { sourceLength, region: region && { start: region.start, end: region.end } } : null;
    });
    return records.length && records.every((record) => record !== null) ? /** @type {LinearDraftRecord[]} */ (records) : [];
  };
  // A Linear card whose record list is still deferred after a Load is read
  // now, through the record selector, as Generate reads it (DESIGN A4).
  /** @param {DrawingState} drawing */
  const awaitLinearRecordDiscovery = async (drawing) => {
    if (drawing.linearSeqs.some((seq) => linearRecordStatusFor(drawing, seq) === 'deferred')) {
      const outcome = refreshLinearRecords(drawing);
      if (outcome instanceof Promise) await outcome;
    }
  };
  // Web-only changes of the shown canvas: a moved diagram, Legend, title, or scale, and Canvas padding.
  const decorationsMoved = () => [state.diagramOffset, state.legendCurrentOffset, state.plotTitleUserOffset, state.lengthBarUserOffset]
    .some((offset) => Boolean(Number(offset.x) || Number(offset.y)));
  /** @param {{ canvasPadding: Record<string, unknown> }} drawing */
  const canvasPadded = (drawing) => (/** @type {const} */ (['top', 'right', 'bottom', 'left'])).some((side) => Number(drawing.canvasPadding[side]));
  /** Whether a drawing shows a Result of its own. @param {{ id: string }} drawing */
  const showsOwnResult = (drawing) => state.generatedDrawingId.value === drawing.id && results.value.length > 0;
  // The committed request of the Result a drawing shows; null while it shows no Result of its own.
  /** @param {{ id: string }} drawing */
  const shownCommittedRequest = (drawing) => (showsOwnResult(drawing) ? getCommittedCanonicalRenderRequest() : null);
  // The record ID the Result a drawing shows names each of its records by (`recordIdsOfResult`); empty while it shows
  // no Result of its own.
  /** @param {{ id: string }} drawing @returns {ReadonlyMap<string, ShownRecordName>} */
  const shownRecordIds = (drawing) => (showsOwnResult(drawing) ? recordIdsOfResult(
    getCommittedCanonicalRenderRequest(), state.trackSlotResolvedGeometry.value, results.value.length
  ) : new Map());
  // Wrap Rows as the shown Result drew it (PD-OI-113), from Python's input
  // summary and the committed request of the drawing whose Result is shown:
  // the field's Auto text and the Auto hints of sizes and ticks read it.
  /** @param {DrawingState} drawing The committed request's stored wrap leaf while this drawing's Result is shown. */
  const shownLinearWrapLeaf = (drawing) => {
    const request = drawing.mode === 'linear' ? shownCommittedRequest(drawing) : null;
    return request?.mode === 'linear' ? requestLinearWrapLeaf(request) : undefined;
  };
  const shownLinearWrapReading = computed(() => {
    const drawing = state.activeDrawing();
    const leaf = shownLinearWrapLeaf(drawing);
    return shownLinearWrap({
      value: drawing.form.linear_wrap_length,
      requested: leaf === undefined ? undefined : wrapFormValue(leaf).value,
      summary: shownInputSummary.value ?? undefined,
      records: linearDraftRecordsOf(drawing)
    });
  });
  return {
    linearRegionNote, linearCardWindowView, linearLayoutCellRows, linearUnitPlacesOf, linearCardRecord, circularRegionNote,
    circularWindowRows, shownInputSummary, drawnLinearRecordLengths, committedLinearRecordsOf, awaitLinearRecordDiscovery, decorationsMoved, canvasPadded, shownCommittedRequest, shownRecordIds,
    shownLinearWrapLeaf, shownLinearWrapReading
  };
};
