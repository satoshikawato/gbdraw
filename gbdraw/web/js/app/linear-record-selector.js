// @ts-check
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { discoveryErrorIsFinal, inputRecordsOperation } from '../services/record-discovery.js';
import { inputSourceFiles, inputTypeSpec } from '../services/input-types.js';
import {
  buildDisambiguatedRecordEntries, formatRecordLength, resolveDisambiguatedRecordSelection
} from '../services/record-options.js';
import { msg, t, tm } from '../utils/i18n.js';
import { errorText } from '../utils/error-text.js';
/** @import { DrawingState, LinearSeq, SessionBusy, WebState } from '../state.js' */
/** @import { InputFileRef } from '../services/history-files.js' */
/** @import { SequenceRecord } from '../services/record-discovery.js' */
/** @import { UserFacingError } from '../utils/error-normalization.js' */

export const AUTOMATIC_RECORD_OPTION_LABEL = msg('Automatic (no explicit selector)');
export const RECORDS_LOADING_LABEL = msg('Loading records...');
// OV-379: a File whose records no read has started for (a loaded Session
// defers discovery) is not inspected; it is not loading.
export const RECORDS_NOT_INSPECTED_LABEL = msg('Not inspected');
export const RECORD_READ_ERROR_LABEL = msg('Records could not be loaded');

/**
 * One record-selector option. `synthetic` marks the entry for a selector the
 * current file does not list.
 * @typedef {{ value: string, label: string, synthetic: boolean }} RecordOption
 */

/** @typedef {'idle' | 'loading' | 'deferred' | 'ready' | 'error'} RecordSelectorStatus */

/**
 * What the selector remembers for one row: the files it read, and how reading went.
 * @typedef {object} RecordSelectorState
 * @property {RecordSelectorStatus} status
 * @property {SequenceRecord[]} records
 * @property {UserFacingError | ''} error The reader's error, or '' when there is none.
 * @property {string} inputType
 * @property {InputFileRef | null} primaryFile
 * @property {InputFileRef | null} pairedFile
 */

/**
 * What `replaceState` takes: it stores a missing field as the empty value, and
 * a falsy error (a thrown null) as ''.
 * @typedef {Omit<Partial<RecordSelectorState>, 'error'> & { status: RecordSelectorStatus, error?: UserFacingError | '' | null }} RecordSelectorStateInput
 */

/** @typedef {Partial<LinearSeq> | null | undefined} LinearSeqLike A row, or a partial one. */
/** @typedef {Pick<DrawingState, 'id' | 'inputType'>} LinearCardDrawing The drawing that holds a card: the drawing object keys the card's state, its input type says which files the card's records come from. */
/** @typedef {Record<string, RecordSelectorState>} RowSelectorStates One drawing's selector states, by row uid. */
/** @typedef {{ suppress?: boolean }} RefreshOptions */
/** @typedef {Pick<DrawingState, 'id' | 'mode' | 'linearSeqs' | 'inputType'>} LinearDiscoveryDrawing The drawing a discovery runs for, over its record rows and input type. */
/** @typedef {SessionBusy | true | undefined} RefreshOutcome What a refresh settles to: a busy outcome, true when it expanded the rows and restarts, or nothing. */

/** @param {LinearSeqLike} seq */
const currentSelectorValue = (seq) => String(seq?.region_record_id ?? '').trim();

export { formatRecordLength };

/**
 * The options of a card's record select. The stored selector names a record as
 * Generate reads it (`resolveDisambiguatedRecordSelection`: by option value,
 * then `#k`, then record ID; OV-3325); that record's option carries the stored
 * value, so the select shows it chosen. A selector that names no single record
 * gets a "not found" entry of its own.
 * @param {readonly SequenceRecord[] | null | undefined} records
 * @param {unknown} [currentValue]
 * @returns {RecordOption[]}
 */
export const buildRecordOptions = (records, currentValue = '') => {
  const normalizedRecords = Array.isArray(records) ? records : [];
  const selectedValue = String(currentValue ?? '').trim();
  const selection = resolveDisambiguatedRecordSelection(normalizedRecords, selectedValue);
  const options = buildDisambiguatedRecordEntries(normalizedRecords).map((record) => {
    const { selector, recordId, usesIndex, value } = record;
    return {
      value: value === selection.record?.value ? selectedValue : value,
      label: `${recordId} (${formatRecordLength(record?.recordLength)})${usesIndex ? ` [${selector}]` : ''}`,
      synthetic: false
    };
  });

  const selectedIsMissing = Boolean(selectedValue) && selection.status !== 'resolved';
  return [
    { value: '', label: tm(AUTOMATIC_RECORD_OPTION_LABEL), synthetic: false },
    ...(selectedIsMissing
      ? [{
          value: selectedValue,
          label: t('{value} (not found in current file)', { value: selectedValue }),
          synthetic: true
        }]
      : []),
    ...options
  ];
};

/** @returns {RecordSelectorState} */
const emptySelectorState = () => ({
  status: 'idle',
  records: [],
  error: '',
  inputType: '',
  primaryFile: null,
  pairedFile: null
});

/** @param {unknown} value */
const sanitizePathSegment = (value) =>
  String(value || 'sequence').replace(/[^A-Za-z0-9_.-]+/g, '_').slice(0, 96) || 'sequence';

/** @param {LinearSeqLike} seq */
const uidFor = (seq) => String(seq?.uid ?? '').trim();

/**
 * Each listed drawing's states of rows it no longer holds go.
 * @param {WeakMap<LinearCardDrawing, RowSelectorStates>} states
 * @param {readonly DrawingState[]} drawings
 */
const purgeRemovedRows = (states, drawings) => drawings.forEach((drawing) => {
  const rows = states.get(drawing) ?? {};
  const activeUids = new Set(drawing.linearSeqs.map(uidFor));
  Object.keys(rows).filter((uid) => !activeUids.has(uid)).forEach((uid) => { delete rows[uid]; });
});

/**
 * A same-ID copy replaced `drawing` (a legacy configuration Load): the copy shows the drawing's record reads; the
 * drawing keeps its own for the Load's rollback.
 * @param {WeakMap<LinearCardDrawing, RowSelectorStates>} states
 * @param {LinearCardDrawing} drawing @param {LinearCardDrawing} copy
 */
const followDrawingCopy = (states, drawing, copy) => {
  const rows = states.get(drawing);
  if (rows) states.set(copy, { ...rows });
};

/**
 * @typedef {Object} LinearRecordReadRequest
 * @property {string} inputType
 * @property {InputFileRef | null} primaryFile
 * @property {InputFileRef | null} pairedFile
 * @property {string} temporaryPathPrefix
 */

/**
 * The members of the `state.js` state that the selector reads; the root passes
 * the whole state. `drawingList()` lists every drawing; a drawing's rows keep
 * their selector state while the drawing holds them; a drawing no longer listed
 * (deleted, or replaced by a Load) keeps its states, and a read that ends after
 * it left lands in its own states but reports no records.
 * @typedef {Pick<WebState, 'drawingList' | 'semanticFileWatchersSuppressed' | 'sessionImportRollbackInProgress'
 *   | 'sessionOperationAvailability'>} LinearRecordSelectorState
 */

/**
 * @typedef {Object} LinearRecordSelectorOptions
 * @property {LinearRecordSelectorState} state
 * @property {(initial: WeakMap<LinearCardDrawing, RowSelectorStates>) => WeakMap<LinearCardDrawing, RowSelectorStates>} reactive
 * @property {(request: LinearRecordReadRequest) => Promise<SequenceRecord[]>} recordReader
 * @property {((discovered: { drawingId: DrawingState['id'], uid: string, records: SequenceRecord[] }) => boolean | void) | null} [onRecordsDiscovered]
 *   The records of one card, for the drawing the discovery ran for. Returns
 *   true when it expanded the rows, which restarts the refresh.
 * @property {{ warn?: (...args: unknown[]) => void }} [logger]
 */

/** @param {LinearRecordSelectorOptions} options */
export const createLinearRecordSelector = ({
  state,
  reactive,
  recordReader,
  onRecordsDiscovered = null,
  logger = console
}) => {
  // OV-2149: each drawing's states, by the drawing object then row uid. A Duplicate
  // keeps its source's row uids, so a uid names a card only inside its drawing.
  // Keyed by the object, a state lives exactly as long as its drawing: a new
  // drawing that takes a deleted drawing's ID starts with none (OV-3892), and
  // the Undo of the Delete brings the drawing back with its own.
  const selectorStateByDrawing = reactive(new WeakMap());
  // Each drawing's discovery: its read generation and the refresh still running.
  /** @type {WeakMap<LinearCardDrawing, number>} */
  const refreshGenerations = new WeakMap();
  /** @type {WeakMap<LinearCardDrawing, { fingerprint: unknown[], promise: Promise<RefreshOutcome> | null }>} */
  const activeRefreshes = new WeakMap();
  const recordDiscoverySuppressed = () => (
    state.semanticFileWatchersSuppressed.value ||
    state.sessionImportRollbackInProgress.value
  );

  /**
   * @param {LinearSeqLike} seq
   * @param {string} inputType the input type of the row's drawing
   */
  const sourceFilesFor = (seq, inputType) => inputSourceFiles(
    inputType, { gb: seq?.gb || null, gff: seq?.gff || null, fasta: seq?.fasta || null }
  );

  /**
   * @param {LinearCardDrawing} drawing the drawing that holds the card
   * @param {LinearSeqLike} seq
   * @returns {RecordSelectorState}
   */
  const stateFor = (drawing, seq) => {
    const uid = uidFor(seq);
    if (!uid) return emptySelectorState();
    // A drawing no refresh has written yet has no entry.
    const stored = selectorStateByDrawing.get(drawing)?.[uid] || emptySelectorState();
    const inputType = drawing.inputType.value;
    const { primaryFile, pairedFile, complete: readable } = sourceFilesFor(seq, inputType);
    if (
      stored.inputType !== inputType ||
      stored.primaryFile !== primaryFile ||
      stored.pairedFile !== pairedFile
    ) {
      // A card this drawing has not read yet (a Duplicate's) shows the settled
      // read of the same files, as its refresh would reuse it.
      const settled = readable ? settledStateFor(inputType, primaryFile, pairedFile) : null;
      if (settled) return settled;
      return {
        ...emptySelectorState(),
        status: readable ? 'deferred' : 'idle',
        inputType,
        primaryFile,
        pairedFile
      };
    }
    return stored;
  };

  /**
   * @param {LinearCardDrawing} drawing
   * @param {string} uid
   * @param {RecordSelectorStateInput} nextState
   */
  const replaceState = (drawing, uid, nextState) => {
    if (!selectorStateByDrawing.has(drawing)) selectorStateByDrawing.set(drawing, {});
    // Set above; read back through the (reactive) map so the write is tracked.
    const rows = /** @type {RowSelectorStates} */ (selectorStateByDrawing.get(drawing));
    rows[uid] = {
      status: nextState.status,
      records: Array.isArray(nextState.records) ? nextState.records : [],
      error: nextState.error || '',
      inputType: String(nextState.inputType || ''),
      primaryFile: nextState.primaryFile || null,
      pairedFile: nextState.pairedFile || null
    };
  };

  // The settled read of these exact source files on any row of any drawing:
  // its records, or an error the reader reported for those bytes. A refresh reuses it.
  /**
   * @param {string} inputType
   * @param {InputFileRef | null} primaryFile
   * @param {InputFileRef | null} pairedFile
   */
  const settledStateFor = (inputType, primaryFile, pairedFile) => state.drawingList()
    .flatMap((drawing) => Object.values(selectorStateByDrawing.get(drawing) ?? {}))
    .find((stored) => (
      stored.inputType === inputType && stored.primaryFile === primaryFile && stored.pairedFile === pairedFile
      && (stored.status === 'ready' || (stored.status === 'error' && discoveryErrorIsFinal(stored.error)))
    )) || null;

  // Each drawing's states of rows it no longer holds go (a drawing's own states go with the drawing).
  const purgeInactiveState = () => purgeRemovedRows(selectorStateByDrawing, state.drawingList());

  /**
   * @param {LinearDiscoveryDrawing} drawing the drawing the discovery runs for
   * @param {{ generation: number, uid: string, primaryFile: InputFileRef | null, pairedFile: InputFileRef | null, inputType: string }} request
   */
  const isCurrentRequest = (drawing, { generation, uid, primaryFile, pairedFile, inputType }) => {
    if (generation !== refreshGenerations.get(drawing)) return false;
    if (recordDiscoverySuppressed()) return false;
    if (drawing.inputType.value !== inputType) return false;
    const currentSeq = drawing.linearSeqs.find((seq) => uidFor(seq) === uid);
    if (!currentSeq) return false;
    const currentFiles = sourceFilesFor(currentSeq, inputType);
    return currentFiles.primaryFile === primaryFile && currentFiles.pairedFile === pairedFile;
  };

  /** @param {LinearDiscoveryDrawing} drawing @param {RefreshOptions} [options] */
  const refreshFingerprint = (drawing, { suppress = false } = {}) => {
    const inputType = drawing.inputType.value;
    return [
      Boolean(suppress || recordDiscoverySuppressed()),
      inputType,
      ...drawing.linearSeqs.flatMap((seq) => {
        const { primaryFile, pairedFile } = sourceFilesFor(seq, inputType);
        return [uidFor(seq), primaryFile, pairedFile];
      })
    ];
  };
  /**
   * @param {readonly unknown[]} left
   * @param {readonly unknown[]} right
   */
  const sameFingerprint = (left, right) => (
    left.length === right.length && left.every((value, index) => Object.is(value, right[index]))
  );

  /**
   * @param {LinearDiscoveryDrawing} drawing the drawing whose cards it reads; its records land only in that drawing
   * @param {RefreshOptions} [options]
   * @returns {Promise<RefreshOutcome>}
   */
  const runRefresh = async (drawing, { suppress = false } = {}) => {
    const generation = (refreshGenerations.get(drawing) ?? 0) + 1;
    refreshGenerations.set(drawing, generation);
    if (suppress || recordDiscoverySuppressed()) return;
    purgeInactiveState();
    if (drawing.mode !== 'linear') return;

    const inputType = drawing.inputType.value;
    /** @type {{ uid: string, primaryFile: InputFileRef, pairedFile: InputFileRef | null }[]} */
    const targets = [];
    for (const seq of drawing.linearSeqs) {
      const uid = uidFor(seq);
      if (!uid) continue;
      const { primaryFile, pairedFile, complete } = sourceFilesFor(seq, inputType);
      if (!primaryFile || !complete) {
        replaceState(drawing, uid, { ...emptySelectorState(), inputType, primaryFile, pairedFile });
        continue;
      }
      if (!settledStateFor(inputType, primaryFile, pairedFile)) {
        replaceState(drawing, uid, {
          status: 'loading', records: [], error: '', inputType, primaryFile, pairedFile
        });
      }
      targets.push({ uid, primaryFile, pairedFile });
    }
    for (const { uid, primaryFile, pairedFile } of targets) {
      const settled = settledStateFor(inputType, primaryFile, pairedFile);
      if (settled?.status === 'error') {
        replaceState(drawing, uid, settled);
        continue;
      }
      try {
        const records = settled ? settled.records : await recordReader({
          inputType,
          primaryFile,
          pairedFile,
          temporaryPathPrefix: `/record-selector-${sanitizePathSegment(uid)}-${generation}`
        });
        if (!isCurrentRequest(drawing, { generation, uid, primaryFile, pairedFile, inputType })) return;
        const busy = state.sessionOperationAvailability();
        if (busy) {
          replaceState(drawing, uid, { status: 'deferred', records: [], error: '', inputType, primaryFile, pairedFile });
          return busy;
        }
        replaceState(drawing, uid, {
          status: 'ready', records, error: '', inputType, primaryFile, pairedFile
        });
        // A drawing no longer listed keeps the result in its own entry (an Undo of its Delete shows it), and
        // its records expand no rows: the discovery finds its drawing by ID, which a new drawing may now hold.
        if (!state.drawingList().some((listed) => listed === drawing)) continue;
        if (onRecordsDiscovered?.({ drawingId: drawing.id, uid, records }) === true) return true;
      } catch (error) {
        if (!isCurrentRequest(drawing, { generation, uid, primaryFile, pairedFile, inputType })) return;
        const busy = state.sessionOperationAvailability();
        if (busy) {
          replaceState(drawing, uid, { status: 'deferred', records: [], error: '', inputType, primaryFile, pairedFile });
          return busy;
        }
        logger.warn?.('Input record discovery failed.', normalizeUserFacingError(error));
        replaceState(drawing, uid, {
          status: 'error',
          records: [],
          error: normalizeUserFacingError(error, { operation: inputRecordsOperation(inputTypeSpec(inputType)), stage: 'helper' }),
          inputType,
          primaryFile,
          pairedFile
        });
      }
    }
  };

  /**
   * The discovery of one drawing's cards (R10: the transition that changed
   * them calls it for that drawing). It shares a refresh still running for
   * the same drawing and sources.
   * @param {LinearDiscoveryDrawing} drawing
   * @param {RefreshOptions} [options]
   * @returns {SessionBusy | Promise<RefreshOutcome>}
   */
  const refresh = (drawing, options = {}) => {
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const fingerprint = refreshFingerprint(drawing, options);
    const activeRefresh = activeRefreshes.get(drawing);
    if (activeRefresh && sameFingerprint(activeRefresh.fingerprint, fingerprint)) {
      // An entry is stored only after its promise is assigned (below).
      return /** @type {Promise<RefreshOutcome>} */ (activeRefresh.promise);
    }
    /** @type {{ fingerprint: unknown[], promise: Promise<RefreshOutcome> | null }} */
    const entry = { fingerprint, promise: null };
    entry.promise = runRefresh(drawing, options).then((expanded) => {
      if (expanded) return refresh(drawing, options);
    }).finally(() => {
      if (activeRefreshes.get(drawing) === entry) activeRefreshes.delete(drawing);
    });
    activeRefreshes.set(drawing, entry);
    // The promise was assigned above.
    return /** @type {Promise<RefreshOutcome>} */ (entry.promise);
  };

  /**
   * @param {LinearCardDrawing} drawing the drawing that holds the card
   * @param {LinearSeqLike} seq
   * @returns {RecordOption[]}
   */
  const optionsFor = (drawing, seq) => {
    const selectorState = stateFor(drawing, seq);
    const currentValue = currentSelectorValue(seq);
    if (selectorState.status === 'ready') {
      return buildRecordOptions(selectorState.records, currentValue);
    }

    const label = selectorState.status === 'loading'
      ? RECORDS_LOADING_LABEL
      : selectorState.status === 'deferred'
        ? RECORDS_NOT_INSPECTED_LABEL
        : selectorState.status === 'error'
          ? RECORD_READ_ERROR_LABEL
          : ''; // No file: the disabled select shows no text (its accessible name says what it is).
    return [{
      value: currentValue,
      label: label ? tm(label) : '',
      synthetic: Boolean(currentValue)
    }];
  };

  // Every reader takes the drawing that holds the card (OV-2149).
  /** @param {LinearCardDrawing} drawing @param {LinearSeqLike} seq */
  const isDisabled = (drawing, seq) => stateFor(drawing, seq).status !== 'ready';
  /** @param {LinearCardDrawing} drawing @param {LinearSeqLike} seq */
  const statusFor = (drawing, seq) => stateFor(drawing, seq).status;
  /** @param {LinearCardDrawing} drawing @param {LinearSeqLike} seq */
  const recordsFor = (drawing, seq) => stateFor(drawing, seq).records.slice();
  /** @param {LinearCardDrawing} drawing @param {LinearSeqLike} seq */
  const errorModelFor = (drawing, seq) => stateFor(drawing, seq).error;
  /** @param {LinearCardDrawing} drawing @param {LinearSeqLike} seq */
  // A stored '' (no error) is no model: it is returned as is.
  const errorFor = (drawing, seq) => {
    const error = errorModelFor(drawing, seq);
    return error && typeof error === 'object' ? errorText(error.summaryParts) : error;
  };
  /**
   * The card's share of its annotation catalog source: the status, records and error model the
   * record select shows, so the catalog says what the card says.
   * @param {LinearCardDrawing} drawing @param {LinearSeqLike} seq
   */
  const catalogStateFor = (drawing, seq) => {
    const { status, records, error } = stateFor(drawing, seq);
    return { status, records: records.slice(), error };
  };
  /** @param {LinearCardDrawing} drawing @param {LinearSeqLike} seq */
  const warningFor = (drawing, seq) => {
    const selectorState = stateFor(drawing, seq);
    const currentValue = currentSelectorValue(seq);
    if (selectorState.status !== 'ready' || !currentValue) return '';
    const available = resolveDisambiguatedRecordSelection(selectorState.records, currentValue).status === 'resolved';
    return available ? '' : t('Selected record was not found in the current file. Choose another record.');
  };

  return {
    selectorStateByDrawing,
    refresh,
    purgeInactiveState,
    followDrawingCopy: (/** @type {LinearCardDrawing} */ drawing, /** @type {LinearCardDrawing} */ copy) => (
      followDrawingCopy(selectorStateByDrawing, drawing, copy)),
    optionsFor,
    statusFor,
    recordsFor,
    isDisabled,
    errorFor,
    errorModelFor,
    catalogStateFor,
    warningFor
  };
};
