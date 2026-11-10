// @ts-check
// The Circular record panel: the record names and the record select's options
// and error, whether a single record's own fields apply, and the Record Order
// rows of a multi-record canvas. The composition root passes the state, the
// mode, the session admission, the shown-drawing readers, and the record
// refresh as ports.
import { circularDiscoveryForInput } from '../services/record-discovery.js';
import { circularDrawingRecordSet, circularRecordEntries, circularRecordSelection, formatRecordLength } from '../services/record-options.js';
import { cardHoldsRecordKey, circularRecordsToDraw, recordKeyExpansion } from '../services/record-draw-selection.js';
import { circularSourceInputs } from '../services/input-types.js';
import { msg, t, tm, tn } from '../utils/i18n.js';
import { circularRecordLayoutOf } from '../services/one-circle-layout.js';
import { RECORDS_NOT_INSPECTED_LABEL } from './linear-record-selector.js';
import { buildDefaultMultiRecordPositions } from './run-analysis.js';
/** @import { DiagramMode, DrawingState, SessionBusy, WebState } from '../state.js' */
/** @import { MultiRecordPosition } from '../services/config.js' */
/** @import { ShownRecordName } from './shown-drawing-readers.js' */

const { computed } = window.Vue;

/**
 * A mode panel's mutating handler: the session's busy answer, nothing when the
 * shown drawing is of another mode, else `action` on the shown drawing.
 * @typedef {<T>(drawingMode: DiagramMode, action: (drawing: DrawingState) => T) => T | SessionBusy | undefined} ShownAction
 */
/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */

/**
 * @typedef {object} CircularRecordPresentationPorts
 * @property {WebState} state
 * @property {{ readonly value: DiagramMode }} mode The shown drawing's mode.
 * @property {{ readonly value: boolean }} semanticMutationAvailable No session operation holds the inputs.
 * @property {() => SessionBusy | null} sessionOperationAvailability
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {() => unknown} refreshCircularRecordOrder The record read of the shown drawing, with its Record Order.
 * @property {(drawing: { id: string }) => ReadonlyMap<string, ShownRecordName>} shownRecordIds The record ID the
 *   Result a drawing shows names each of its records by, by record key; empty while it shows no Result of its own
 *   (app/shown-drawing-readers.js).
 */

// The file name a record's label begins with: its card's file when the drawing has several source cards.
/** @param {DrawingState} drawing @returns {(uid: string) => string} */
const cardFileNames = (drawing) => {
  const cards = circularSourceInputs(drawing.inputType.value, drawing.circularSources);
  const names = new Map(cards.length > 1 ? cards.map((card) => [card.uid, String(card.primaryFile?.name || '')]) : []);
  return (uid) => names.get(uid) || '';
};

// A record's name, read or not (D-98, D-100): its card's file name when the drawing
// has several cards (`cardFileNames`), its record ID, and its `#k` when that ID
// repeats in its file or is unsafe (`usesIndex`); the Record select writes the
// record's length between the ID and the `#k`.
/** @param {string} fileName @param {string} recordId */
const recordName = (fileName, recordId) => `${fileName ? `${fileName}: ` : ''}${recordId}`;
/** @param {{ usesIndex: boolean, selector: string }} entry */
const recordIndexSuffix = (entry) => (entry.usesIndex ? ` [${entry.selector}]` : '');

/**
 * The name of a value no entry names (`recordEntryOf`): a record identity by its
 * card's file name, with `#k` for the k-th record of a file; any other value (a
 * kept `#N`, a record ID) as it is.
 * @param {DrawingState} drawing
 * @param {string} value
 */
const unnamedRecordName = (drawing, value) => {
  const card = drawing.circularSources.find((source) => cardHoldsRecordKey(source.uid, value));
  if (!card) return value;
  const [input] = circularSourceInputs(drawing.inputType.value, [card]);
  const fileName = String(input.primaryFile?.name || '') || card.uid;
  const k = recordKeyExpansion(value, card.uid);
  return k === null ? fileName : `${fileName} #${k}`;
};

/**
 * The record select, its inspection, and the single-record fields of the shown Circular drawing.
 * @param {CircularRecordPresentationPorts} ports
 */
export const createCircularRecordPresentation = ({
  state, mode, semanticMutationAvailable, sessionOperationAvailability, shownDrawingOfMode, refreshCircularRecordOrder,
  shownRecordIds
}) => {
  /**
   * The named entry of a record identity (`<uid>`, `<uid>:<k>`), read or not; the
   * one owner of the parts a Record select or Record Order label names a record
   * by. A read record's entry is its `circularRecordEntries` entry. An unread one
   * (OV-2949, OV-2959) is named from the records the shown Result names
   * (`shownRecordIds`) the same way; as the Result names only the records it drew,
   * the k-th record of a file keeps its `#k` unless its ID suffices
   * (`ShownRecordName`): the committed request selected it by its record ID, which
   * the Web writes only for an ID that is unique and safe in its file, or drew its
   * whole file, so a repeat shows among the named records (reviews 1 and 2 of
   * #1778). Null for a record the Result does not
   * name and for a value that is no identity.
   * @param {DrawingState} drawing
   * @param {string} value
   */
  const recordEntryOf = (drawing, value) => {
    const read = circularRecordEntries(drawing.circularRecordList.value).find((record) => record.value === value);
    if (read) return read;
    const card = drawing.circularSources.find((source) => cardHoldsRecordKey(source.uid, value));
    if (!card) return null;
    const named = [...shownRecordIds(drawing)].filter(([key]) => cardHoldsRecordKey(card.uid, key))
      .map(([identity, { recordId, idSuffices }]) => {
        const k = recordKeyExpansion(identity, card.uid);
        return { uid: card.uid, identity, k: k ?? 1, selector: `#${k ?? 1}`, record_id: recordId, provenUnique: k === null || idSuffices };
      });
    const entry = circularRecordEntries(named).find((record) => record.value === value);
    return entry ? { ...entry, usesIndex: entry.usesIndex || !entry.provenUnique } : null;
  };
  // The Record select's name of a value whose record it does not list (its synthetic option and its error).
  /** @param {DrawingState} drawing @param {string} value */
  const selectedRecordName = (drawing, value) => {
    const entry = recordEntryOf(drawing, value);
    return entry ? `${recordName(cardFileNames(drawing)(String(entry.uid)), entry.recordId)}${recordIndexSuffix(entry)}`
      : unnamedRecordName(drawing, value);
  };
  // A Record Order row names its record by identity, read or not, from its entry: with
  // several cards as the Record select does, with one card as `#k (<record ID>)`.
  /** @param {unknown} identity */
  const getCircularRecordOrderLabel = (identity) => {
    const normalized = String(identity || '').trim();
    if (!normalized) return '';
    const drawing = state.activeDrawing();
    const entry = recordEntryOf(drawing, normalized);
    if (!entry) return unnamedRecordName(drawing, normalized);
    const fileName = cardFileNames(drawing)(String(entry.uid));
    return fileName ? `${recordName(fileName, entry.recordId)}${recordIndexSuffix(entry)}` : `${entry.selector} (${entry.recordId})`;
  };

  const circularRecordDiscoveryState = computed(() => circularDiscoveryForInput(state.activeDrawing()));
  const circularRecordInspectionEnabled = computed(() => semanticMutationAvailable.value
    && mode.value === 'circular'
    && circularRecordDiscoveryState.value.hasInput
    && circularRecordDiscoveryState.value.status !== 'loading'
    && !state.semanticFileWatchersSuppressed.value);
  const inspectCircularSourceRecords = () => sessionOperationAvailability() || (circularRecordInspectionEnabled.value
    ? refreshCircularRecordOrder()
    : Promise.resolve({ status: 'unavailable', reason: msg('Source inspection is unavailable.') }));
  const circularRecordPresentationEntries = () => circularRecordEntries(circularRecordDiscoveryState.value.records);

  // OV-3340: the records the next Generate draws (`circularDrawingRecordSet`, the request's own owner): the record
  // the Record select names alone, else the ON records.
  const circularDrawnRecordCount = computed(() => {
    const drawing = shownDrawingOfMode('circular');
    return drawing ? circularDrawingRecordSet(circularRecordDiscoveryState.value.records, drawing).records.length : 0;
  });
  const circularRecordPresentationOptions = computed(() => {
    const drawing = state.activeDrawing();
    const entries = circularRecordPresentationEntries();
    const fileNameOf = cardFileNames(drawing);
    const current = String(drawing.form.circular_record_selector || '').trim();
    const selection = circularRecordSelection(entries, current);
    const shown = selectedRecordName(drawing, current);
    // Without an inspected catalog the saved selector is unverified, not missing.
    const inspected = circularRecordDiscoveryState.value.status === 'ready';
    // The automatic option describes the all-records choice: the ON records it draws.
    const drawnCount = circularRecordsToDraw(entries, drawing.recordsOff).length;
    const layout = circularRecordLayoutOf(drawing);
    // OV-379: a source whose records are not read says so, as the Linear select does.
    const automaticLabel = circularRecordDiscoveryState.value.status === 'deferred'
      ? tm(RECORDS_NOT_INSPECTED_LABEL)
      // OV-2839, D-11: the canvas and one circle draw every ON record in one figure; the label names the layout.
      : layout === 'one_circle' && entries.length > 1
        ? (drawnCount < entries.length ? tn(entries.length, '{drawn} of {count} record (One circle)', '{drawn} of {count} records (One circle)', { drawn: drawnCount }) : t('All records (One circle)'))
      : layout === 'grid' && entries.length > 1
        ? (drawnCount < entries.length ? tn(entries.length, '{drawn} of {count} record (Multi-Record Canvas)', '{drawn} of {count} records (Multi-Record Canvas)', { drawn: drawnCount }) : t('All records (Multi-Record Canvas)'))
        : entries.length > 1 || drawing.adv.circular_grouping_intent === 'batch'
          ? (drawnCount < entries.length ? tn(entries.length, '{drawn} of {count} record (separate diagrams)', '{drawn} of {count} records (separate diagrams)', { drawn: drawnCount }) : t('All records (separate diagrams)'))
          : t('Automatic (single record)');
    return [
      { value: '', label: automaticLabel, synthetic: false },
      ...(current && selection.status !== 'resolved'
        ? [{
            value: current,
            label: !inspected ? t('{selector} (not inspected)', { selector: shown })
              : selection.status === 'ambiguous' ? t('{selector} (ambiguous)', { selector: shown })
                : t('{selector} (not found)', { selector: shown }),
            synthetic: true
          }]
        : []),
      ...entries.map((record) => ({
        value: record.value,
        label: `${recordName(fileNameOf(String(record.uid)), record.recordId)} (${formatRecordLength(record.recordLength)})${recordIndexSuffix(record)}`,
        synthetic: false
      }))
    ];
  });

  const circularRecordPresentationError = computed(() => {
    const drawing = state.activeDrawing();
    const current = String(drawing.form.circular_record_selector || '').trim();
    if (!current || circularRecordDiscoveryState.value.status !== 'ready') return '';
    const selection = circularRecordSelection(circularRecordPresentationEntries(), current);
    if (selection.status === 'ambiguous') {
      return t('Record selector \'{selector}\' is ambiguous in the current input.', { selector: selectedRecordName(drawing, current) });
    }
    if (selection.status === 'missing') {
      return t('No record in the current input matches the record selector \'{selector}\'.', { selector: selectedRecordName(drawing, current) });
    }
    return '';
  });

  // Why a Circular record's own fields (crop, orientation, title lines) are
  // unavailable, in the order the sidebar explains it; '' when they are not.
  const circularOwnRecordBlock = computed(() => {
    const drawing = state.activeDrawing();
    const layout = circularRecordLayoutOf(drawing);
    if (layout !== 'separate') return layout === 'one_circle' ? 'one-circle' : 'multi-record-canvas';
    if (circularRecordDiscoveryState.value.status !== 'ready') return 'not-inspected';
    const entries = circularRecordPresentationEntries();
    const selection = circularRecordSelection(entries, drawing.form.circular_record_selector);
    return drawing.adv.circular_grouping_intent !== 'batch'
      && (selection.status === 'resolved' || (selection.status === 'unspecified' && entries.length === 1))
      ? '' : 'no-single-record';
  });
  const circularSingleRecordPresentationEnabled = computed(() => mode.value === 'circular' && !circularOwnRecordBlock.value);
  // The shown drawing's record layout (`circularRecordLayoutOf`), which the layout checkboxes and texts name.
  const circularRecordLayout = computed(() => circularRecordLayoutOf(state.activeDrawing()));
  return {
    getCircularRecordOrderLabel, circularRecordDiscoveryState, circularRecordInspectionEnabled, inspectCircularSourceRecords,
    circularDrawnRecordCount, circularRecordPresentationOptions, circularRecordPresentationError, circularOwnRecordBlock,
    circularSingleRecordPresentationEnabled, circularRecordLayout
  };
};

/** @param {DrawingState} drawing */
export const buildDefaultCircularRecordPositions = (drawing) => buildDefaultMultiRecordPositions(
  drawing.circularRecordList.value.map((entry) => entry.identity)
);

/**
 * @param {DrawingState} drawing
 * @param {MultiRecordPosition | null | undefined} position
 */
const circularRecordRowOf = (drawing, position) => {
  const rowValue = Number(position?.row);
  const maxRow = Math.max(1, Array.isArray(drawing.adv.multi_record_positions) ? drawing.adv.multi_record_positions.length : 1);
  if (!Number.isInteger(rowValue) || rowValue <= 0) return 1;
  return Math.min(rowValue, maxRow);
};

/** @param {DrawingState} drawing */
const sortCircularRecordPositionsByRow = (drawing) => {
  if (!Array.isArray(drawing.adv.multi_record_positions) || drawing.adv.multi_record_positions.length <= 1) return;
  const sorted = drawing.adv.multi_record_positions
    .map((/** @type {MultiRecordPosition} */ entry, /** @type {number} */ index) => ({ ...entry, __index: index }))
    .sort((/** @type {MultiRecordPosition & { __index: number }} */ left, /** @type {MultiRecordPosition & { __index: number }} */ right) => {
      const leftRow = Number(left.row);
      const rightRow = Number(right.row);
      if (leftRow !== rightRow) return leftRow - rightRow;
      return left.__index - right.__index;
    })
    .map((/** @type {MultiRecordPosition & { __index: number }} */ { __index, ...entry }) => entry);
  drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length, ...sorted);
};

/** @param {DrawingState} drawing @param {number} index */
const canMoveCircularRecordOrderUpIn = (drawing, index) => {
  const idx = Number(index);
  if (!Number.isInteger(idx) || idx <= 0 || idx >= drawing.adv.multi_record_positions.length) return false;
  const currentRow = circularRecordRowOf(drawing, drawing.adv.multi_record_positions[idx]);
  const prevRow = circularRecordRowOf(drawing, drawing.adv.multi_record_positions[idx - 1]);
  return currentRow === prevRow;
};

/** @param {DrawingState} drawing @param {number} index */
const canMoveCircularRecordOrderDownIn = (drawing, index) => {
  const idx = Number(index);
  if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.adv.multi_record_positions.length - 1) return false;
  const currentRow = circularRecordRowOf(drawing, drawing.adv.multi_record_positions[idx]);
  const nextRow = circularRecordRowOf(drawing, drawing.adv.multi_record_positions[idx + 1]);
  return currentRow === nextRow;
};

/**
 * @typedef {object} CircularRecordOrderPorts
 * @property {ShownAction} shownAction
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 */

/**
 * The Record Order rows of the shown Circular drawing: each record's row, and its moves within the row.
 * @param {CircularRecordOrderPorts} ports
 */
export const createCircularRecordOrder = ({ shownAction, shownDrawingOfMode }) => {
  /** @param {MultiRecordPosition | null | undefined} position */
  const getCircularRecordRow = (position) => {
    const drawing = shownDrawingOfMode('circular');
    return drawing ? circularRecordRowOf(drawing, position) : 1;
  };

  const getCircularRecordRowOptions = () => {
    const drawing = shownDrawingOfMode('circular');
    if (!drawing) return [1];
    const count = Array.isArray(drawing.adv.multi_record_positions) ? drawing.adv.multi_record_positions.length : 0;
    const maxRow = Math.max(1, count);
    return Array.from({ length: maxRow }, (_unused, index) => index + 1);
  };

  /** @param {number} index @param {unknown} rowValue */
  const setCircularRecordRow = (index, rowValue) => shownAction('circular', (drawing) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.adv.multi_record_positions.length) return;
    const target = drawing.adv.multi_record_positions[idx];
    if (!target || typeof target !== 'object' || Array.isArray(target)) return;
    const maxRow = Math.max(1, drawing.adv.multi_record_positions.length);
    const parsedRow = Number(rowValue);
    const normalizedRow = Number.isInteger(parsedRow) && parsedRow > 0 ? Math.min(parsedRow, maxRow) : 1;
    target.row = normalizedRow;
    sortCircularRecordPositionsByRow(drawing);
  });

  /** @param {number} index */
  const canMoveCircularRecordOrderUp = (index) => {
    const drawing = shownDrawingOfMode('circular');
    return drawing ? canMoveCircularRecordOrderUpIn(drawing, index) : false;
  };

  /** @param {number} index */
  const canMoveCircularRecordOrderDown = (index) => {
    const drawing = shownDrawingOfMode('circular');
    return drawing ? canMoveCircularRecordOrderDownIn(drawing, index) : false;
  };

  /** @param {number} index */
  const moveCircularRecordOrderUp = (index) => shownAction('circular', (drawing) => {
    const idx = Number(index);
    if (!canMoveCircularRecordOrderUpIn(drawing, idx)) return;
    const next = [...drawing.adv.multi_record_positions];
    const temp = next[idx - 1];
    next[idx - 1] = next[idx];
    next[idx] = temp;
    drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length, ...next);
  });

  /** @param {number} index */
  const moveCircularRecordOrderDown = (index) => shownAction('circular', (drawing) => {
    const idx = Number(index);
    if (!canMoveCircularRecordOrderDownIn(drawing, idx)) return;
    const next = [...drawing.adv.multi_record_positions];
    const temp = next[idx + 1];
    next[idx + 1] = next[idx];
    next[idx] = temp;
    drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length, ...next);
  });

  const resetCircularRecordOrder = () => shownAction('circular', (drawing) => {
    const defaults = buildDefaultCircularRecordPositions(drawing);
    drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length, ...defaults);
  });
  return {
    getCircularRecordRow, getCircularRecordRowOptions, setCircularRecordRow, canMoveCircularRecordOrderUp,
    canMoveCircularRecordOrderDown, moveCircularRecordOrderUp, moveCircularRecordOrderDown, resetCircularRecordOrder
  };
};
