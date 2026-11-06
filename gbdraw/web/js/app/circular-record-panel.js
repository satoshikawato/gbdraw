// @ts-check
// The Circular record panel: the record names and the record select's options
// and error, whether a single record's own fields apply, and the Record Order
// rows of a multi-record canvas. The composition root passes the state, the
// mode, the session admission, the shown-drawing readers, and the record
// refresh as ports.
import { circularDiscoveryForInput } from '../services/record-discovery.js';
import {
  buildDisambiguatedRecordEntries,
  formatRecordLength,
  resolveDisambiguatedRecordSelection
} from '../services/record-options.js';
import { circularRecordsToDraw } from '../services/record-draw-selection.js';
import { msg, t, tm, tn } from '../utils/i18n.js';
import { RECORDS_NOT_INSPECTED_LABEL } from './linear-record-selector.js';
import { buildDefaultMultiRecordPositions } from './run-analysis.js';
/** @import { DiagramMode, DrawingState, SessionBusy, WebState } from '../state.js' */
/** @import { SequenceRecordEntry } from '../services/record-discovery.js' */
/** @import { MultiRecordPosition } from '../services/config.js' */

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
 */

/**
 * The record select, its inspection, and the single-record fields of the shown Circular drawing.
 * @param {CircularRecordPresentationPorts} ports
 */
export const createCircularRecordPresentation = ({
  state, mode, semanticMutationAvailable, sessionOperationAvailability, shownDrawingOfMode, refreshCircularRecordOrder
}) => {
  /** @param {unknown} selector */
  const getCircularRecordOrderLabel = (selector) => {
    const normalized = String(selector || '').trim();
    if (!normalized) return '';
    const records = state.activeDrawing().circularRecordList.value;
    const matched = records.find((/** @type {{ selector?: unknown }} */ entry) => String(entry?.selector || '').trim() === normalized);
    if (!matched) return normalized;
    return `${normalized} (${String(matched.record_id || '').trim() || t('Unknown')})`;
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
  const circularRecordPresentationEntries = () => buildDisambiguatedRecordEntries(
    circularRecordDiscoveryState.value.records.map(
      (/** @type {SequenceRecordEntry} */ record) => ({
        ...record,
        recordId: record?.record_id ?? record?.recordId,
        recordLength: record?.record_length ?? record?.recordLength
      })
    )
  );

  const circularDrawnRecordCount = computed(() => {
    const drawing = shownDrawingOfMode('circular');
    return drawing ? circularRecordsToDraw(circularRecordPresentationEntries(), drawing.recordsOff).length : 0;
  });
  const circularRecordPresentationOptions = computed(() => {
    const drawing = state.activeDrawing();
    const entries = circularRecordPresentationEntries();
    const current = String(drawing.form.circular_record_selector || '').trim();
    const selection = resolveDisambiguatedRecordSelection(entries, current);
    // Without an inspected catalog the saved selector is unverified, not missing.
    const inspected = circularRecordDiscoveryState.value.status === 'ready';
    const drawnCount = circularDrawnRecordCount.value;
    // OV-379: a source whose records are not read says so, as the Linear select does.
    const automaticLabel = circularRecordDiscoveryState.value.status === 'deferred'
      ? tm(RECORDS_NOT_INSPECTED_LABEL)
      : entries.length > 1 || drawing.adv.circular_grouping_intent === 'batch'
        ? (drawnCount < entries.length ? tn(entries.length, '{drawn} of {count} record (separate diagrams)', '{drawn} of {count} records (separate diagrams)', { drawn: drawnCount }) : t('All records (separate diagrams)'))
        : t('Automatic (single record)');
    return [
      { value: '', label: automaticLabel, synthetic: false },
      ...(current && selection.status !== 'resolved'
        ? [{
            value: current,
            label: !inspected ? t('{selector} (not inspected)', { selector: current })
              : selection.status === 'ambiguous' ? t('{selector} (ambiguous)', { selector: current })
                : t('{selector} (not found)', { selector: current }),
            synthetic: true
          }]
        : []),
      ...entries.map((record) => ({
        value: record.value,
        label: `${record.recordId} (${formatRecordLength(record.recordLength)})${record.usesIndex ? ` [${record.selector}]` : ''}`,
        synthetic: false
      }))
    ];
  });

  const circularRecordPresentationError = computed(() => {
    const drawing = state.activeDrawing();
    const current = String(drawing.form.circular_record_selector || '').trim();
    if (!current || circularRecordDiscoveryState.value.status !== 'ready') return '';
    const selection = resolveDisambiguatedRecordSelection(
      circularRecordPresentationEntries(),
      current
    );
    if (selection.status === 'ambiguous') {
      return t('Record selector \'{selector}\' is ambiguous in the current input.', { selector: current });
    }
    if (selection.status === 'missing') {
      return t('No record in the current input matches the record selector \'{selector}\'.', { selector: current });
    }
    return '';
  });

  // Why a Circular record's own fields (crop, orientation, title lines) are
  // unavailable, in the order the sidebar explains it; '' when they are not.
  const circularOwnRecordBlock = computed(() => {
    const drawing = state.activeDrawing();
    if (drawing.form.multi_record_canvas) return 'multi-record-canvas';
    if (circularRecordDiscoveryState.value.status !== 'ready') return 'not-inspected';
    const entries = circularRecordPresentationEntries();
    const selection = resolveDisambiguatedRecordSelection(entries, drawing.form.circular_record_selector);
    return drawing.adv.circular_grouping_intent !== 'batch'
      && (selection.status === 'resolved' || (selection.status === 'unspecified' && entries.length === 1))
      ? '' : 'no-single-record';
  });
  const circularSingleRecordPresentationEnabled = computed(() => mode.value === 'circular' && !circularOwnRecordBlock.value);
  return {
    getCircularRecordOrderLabel, circularRecordDiscoveryState, circularRecordInspectionEnabled, inspectCircularSourceRecords,
    circularDrawnRecordCount, circularRecordPresentationOptions, circularRecordPresentationError, circularOwnRecordBlock,
    circularSingleRecordPresentationEnabled
  };
};

/** @param {DrawingState} drawing */
export const buildDefaultCircularRecordPositions = (drawing) => buildDefaultMultiRecordPositions(
  drawing.circularRecordList.value
    .map((/** @type {{ selector?: unknown }} */ entry) => entry?.selector)
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
