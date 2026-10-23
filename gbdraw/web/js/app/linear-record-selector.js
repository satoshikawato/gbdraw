// @ts-check
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { discoveryErrorIsFinal } from '../services/record-discovery.js';
import { buildDisambiguatedRecordEntries, formatRecordLength } from '../services/record-options.js';
import { msg, t, tm } from '../utils/i18n.js';
import { errorText } from '../utils/error-text.js';
/** @import { DrawingState, LinearSeq, DiagramMode, SessionBusy } from '../state.js' */
/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { InputFileRef } from '../services/history-files.js' */
/** @import { SequenceRecord } from '../services/record-discovery.js' */
/** @import { UserFacingError } from '../utils/error-normalization.js' */

export const AUTOMATIC_RECORD_OPTION_LABEL = msg('Automatic (no explicit selector)');
export const RECORDS_LOADING_LABEL = msg('Loading records...');
// OV-379: a File whose records no read has started for (a loaded Session
// defers discovery) is not inspected; it is not loading.
export const RECORDS_NOT_INSPECTED_LABEL = msg('Not inspected');
export const RECORD_FILE_REQUIRED_LABEL = msg('Upload a sequence file first');
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
/** @typedef {{ suppress?: boolean }} RefreshOptions */
/** @typedef {SessionBusy | true | undefined} RefreshOutcome What a refresh settles to: a busy outcome, true when it expanded the rows and restarts, or nothing. */

/** @param {LinearSeqLike} seq */
const currentSelectorValue = (seq) => String(seq?.region_record_id ?? '').trim();

export { formatRecordLength };

/**
 * @param {readonly SequenceRecord[] | null | undefined} records
 * @param {unknown} [currentValue]
 * @returns {RecordOption[]}
 */
export const buildRecordOptions = (records, currentValue = '') => {
  const normalizedRecords = Array.isArray(records) ? records : [];
  const options = buildDisambiguatedRecordEntries(normalizedRecords).map((record) => {
    const { selector, recordId, usesIndex, value } = record;
    return {
      value,
      label: `${recordId} (${formatRecordLength(record?.recordLength)})${usesIndex ? ` [${selector}]` : ''}`,
      synthetic: false
    };
  });

  const selectedValue = String(currentValue ?? '').trim();
  const selectedIsMissing = Boolean(selectedValue) && !options.some((option) => option.value === selectedValue);
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

/**
 * @typedef {Object} LinearRecordReadRequest
 * @property {string} inputType
 * @property {InputFileRef | null} primaryFile
 * @property {InputFileRef | null} pairedFile
 * @property {string} temporaryPathPrefix
 */

/**
 * The members of the `state.js` state that the selector reads; the root passes
 * the whole state. The two suppression flags and `sessionOperationAvailability`
 * are absent in a test state.
 * @typedef {object} LinearRecordSelectorState
 * @property {{ value: string }} lInputType
 * @property {{ value: DiagramMode }} mode
 * @property {() => Pick<DrawingState, 'id' | 'mode'>} activeDrawing the shown drawing, which a discovery runs for
 * @property {LinearSeq[]} linearSeqs
 * @property {{ value: boolean }} [semanticFileWatchersSuppressed]
 * @property {{ value: boolean }} [sessionImportRollbackInProgress]
 * @property {() => SessionBusy | null} [sessionOperationAvailability]
 */

/**
 * @typedef {Object} LinearRecordSelectorOptions
 * @property {LinearRecordSelectorState} state
 * @property {(initial: Record<string, RecordSelectorState>) => Record<string, RecordSelectorState>} reactive
 * @property {(request: LinearRecordReadRequest) => Promise<SequenceRecord[]>} recordReader
 * @property {((discovered: { drawingId: DrawingId, uid: string, records: SequenceRecord[] }) => boolean | void) | null} [onRecordsDiscovered]
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
  const selectorStateByUid = reactive({});
  let refreshGeneration = 0;
  /** @type {{ fingerprint: unknown[], promise: Promise<RefreshOutcome> | null } | null} */
  let activeRefresh = null;
  const recordDiscoverySuppressed = () => Boolean(
    state.semanticFileWatchersSuppressed?.value ||
    state.sessionImportRollbackInProgress?.value
  );

  /** @param {LinearSeqLike} seq */
  const uidFor = (seq) => String(seq?.uid ?? '').trim();
  /**
   * @param {LinearSeqLike} seq
   * @param {string} [inputType]
   * @returns {{ primaryFile: InputFileRef | null, pairedFile: InputFileRef | null }}
   */
  const sourceFilesFor = (seq, inputType = state.lInputType.value) => (
    inputType === 'gff'
      ? { primaryFile: seq?.gff || null, pairedFile: seq?.fasta || null }
      : { primaryFile: seq?.gb || null, pairedFile: null }
  );

  /**
   * @param {LinearSeqLike} seq
   * @returns {RecordSelectorState}
   */
  const stateFor = (seq) => {
    const uid = uidFor(seq);
    const stored = (uid && selectorStateByUid[uid]) || emptySelectorState();
    const inputType = state.lInputType.value;
    const { primaryFile, pairedFile } = sourceFilesFor(seq, inputType);
    if (
      stored.inputType !== inputType ||
      stored.primaryFile !== primaryFile ||
      stored.pairedFile !== pairedFile
    ) {
      return {
        ...emptySelectorState(),
        status: primaryFile && (inputType !== 'gff' || pairedFile) ? 'deferred' : 'idle',
        inputType,
        primaryFile,
        pairedFile
      };
    }
    return stored;
  };

  /**
   * @param {string} uid
   * @param {RecordSelectorStateInput} nextState
   */
  const replaceState = (uid, nextState) => {
    selectorStateByUid[uid] = {
      status: nextState.status,
      records: Array.isArray(nextState.records) ? nextState.records : [],
      error: nextState.error || '',
      inputType: String(nextState.inputType || ''),
      primaryFile: nextState.primaryFile || null,
      pairedFile: nextState.pairedFile || null
    };
  };

  // The settled read of these exact source files on any row: its records, or
  // an error the reader reported for those bytes. A refresh reuses it.
  /**
   * @param {string} inputType
   * @param {InputFileRef | null} primaryFile
   * @param {InputFileRef | null} pairedFile
   */
  const settledStateFor = (inputType, primaryFile, pairedFile) => Object.values(selectorStateByUid)
    .find((stored) => (
      stored.inputType === inputType && stored.primaryFile === primaryFile && stored.pairedFile === pairedFile
      && (stored.status === 'ready' || (stored.status === 'error' && discoveryErrorIsFinal(stored.error)))
    )) || null;

  const purgeInactiveState = () => {
    const activeUids = new Set(state.linearSeqs.map(uidFor).filter(Boolean));
    Object.keys(selectorStateByUid).forEach((uid) => {
      if (!activeUids.has(uid)) delete selectorStateByUid[uid];
    });
  };

  /** @param {{ generation: number, uid: string, primaryFile: InputFileRef | null, pairedFile: InputFileRef | null, inputType: string }} request */
  const isCurrentRequest = ({ generation, uid, primaryFile, pairedFile, inputType }) => {
    if (generation !== refreshGeneration) return false;
    if (recordDiscoverySuppressed()) return false;
    if (state.mode.value !== 'linear' || state.lInputType.value !== inputType) return false;
    const currentSeq = state.linearSeqs.find((seq) => uidFor(seq) === uid);
    if (!currentSeq) return false;
    const currentFiles = sourceFilesFor(currentSeq, inputType);
    return currentFiles.primaryFile === primaryFile && currentFiles.pairedFile === pairedFile;
  };

  /** @param {RefreshOptions} [options] */
  const refreshFingerprint = ({ suppress = false } = {}) => [
    Boolean(suppress || recordDiscoverySuppressed()),
    state.mode.value,
    state.lInputType.value,
    ...state.linearSeqs.flatMap((seq) => {
      const { primaryFile, pairedFile } = sourceFilesFor(seq);
      return [uidFor(seq), primaryFile, pairedFile];
    })
  ];
  /**
   * @param {readonly unknown[]} left
   * @param {readonly unknown[]} right
   */
  const sameFingerprint = (left, right) => (
    left.length === right.length && left.every((value, index) => Object.is(value, right[index]))
  );

  /**
   * @param {RefreshOptions} [options]
   * @returns {Promise<RefreshOutcome>}
   */
  const runRefresh = async ({ suppress = false } = {}) => {
    const generation = ++refreshGeneration;
    if (suppress || recordDiscoverySuppressed()) return;
    purgeInactiveState();

    // An inactive Linear source keeps its records, as an inactive Circular one
    // does. The drawing this discovery runs for is pinned before its first
    // read: its records land in that drawing.
    const drawing = state.activeDrawing();
    if (drawing.mode !== 'linear') return;

    const inputType = state.lInputType.value;
    /** @type {{ uid: string, primaryFile: InputFileRef, pairedFile: InputFileRef | null }[]} */
    const targets = [];
    for (const seq of state.linearSeqs) {
      const uid = uidFor(seq);
      if (!uid) continue;
      const { primaryFile, pairedFile } = sourceFilesFor(seq, inputType);
      if (!primaryFile || (inputType === 'gff' && !pairedFile)) {
        replaceState(uid, { ...emptySelectorState(), inputType, primaryFile, pairedFile });
        continue;
      }
      if (!settledStateFor(inputType, primaryFile, pairedFile)) {
        replaceState(uid, {
          status: 'loading', records: [], error: '', inputType, primaryFile, pairedFile
        });
      }
      targets.push({ uid, primaryFile, pairedFile });
    }
    for (const { uid, primaryFile, pairedFile } of targets) {
      const settled = settledStateFor(inputType, primaryFile, pairedFile);
      if (settled?.status === 'error') {
        replaceState(uid, settled);
        continue;
      }
      try {
        const records = settled ? settled.records : await recordReader({
          inputType,
          primaryFile,
          pairedFile,
          temporaryPathPrefix: `/record-selector-${sanitizePathSegment(uid)}-${generation}`
        });
        if (!isCurrentRequest({ generation, uid, primaryFile, pairedFile, inputType })) return;
        const busy = state.sessionOperationAvailability?.();
        if (busy) {
          replaceState(uid, { status: 'deferred', records: [], error: '', inputType, primaryFile, pairedFile });
          return busy;
        }
        replaceState(uid, {
          status: 'ready', records, error: '', inputType, primaryFile, pairedFile
        });
        if (onRecordsDiscovered?.({ drawingId: drawing.id, uid, records }) === true) return true;
      } catch (error) {
        if (!isCurrentRequest({ generation, uid, primaryFile, pairedFile, inputType })) return;
        const busy = state.sessionOperationAvailability?.();
        if (busy) {
          replaceState(uid, { status: 'deferred', records: [], error: '', inputType, primaryFile, pairedFile });
          return busy;
        }
        logger.warn?.('Input record discovery failed.', normalizeUserFacingError(error));
        replaceState(uid, {
          status: 'error',
          records: [],
          error: normalizeUserFacingError(error, { operation: inputType === 'gff' ? 'listGffFastaRecords' : 'listSequenceRecords', stage: 'helper' }),
          inputType,
          primaryFile,
          pairedFile
        });
      }
    }
  };

  /**
   * @param {RefreshOptions} [options]
   * @returns {SessionBusy | Promise<RefreshOutcome>}
   */
  const refresh = (options = {}) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const fingerprint = refreshFingerprint(options);
    if (activeRefresh && sameFingerprint(activeRefresh.fingerprint, fingerprint)) {
      // An entry is stored only after its promise is assigned (below).
      return /** @type {Promise<RefreshOutcome>} */ (activeRefresh.promise);
    }
    /** @type {{ fingerprint: unknown[], promise: Promise<RefreshOutcome> | null }} */
    const entry = { fingerprint, promise: null };
    entry.promise = runRefresh(options).then((expanded) => {
      if (expanded) return refresh(options);
    }).finally(() => {
      if (activeRefresh === entry) activeRefresh = null;
    });
    activeRefresh = entry;
    // The promise was assigned above.
    return /** @type {Promise<RefreshOutcome>} */ (entry.promise);
  };

  /**
   * @param {LinearSeqLike} seq
   * @returns {RecordOption[]}
   */
  const optionsFor = (seq) => {
    const selectorState = stateFor(seq);
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
          : RECORD_FILE_REQUIRED_LABEL;
    return [{
      value: currentValue,
      label: tm(label),
      synthetic: Boolean(currentValue)
    }];
  };

  /** @param {LinearSeqLike} seq */
  const isDisabled = (seq) => stateFor(seq).status !== 'ready';
  /** @param {LinearSeqLike} seq */
  const statusFor = (seq) => stateFor(seq).status;
  /** @param {LinearSeqLike} seq */
  const recordsFor = (seq) => stateFor(seq).records.slice();
  /** @param {LinearSeqLike} seq */
  const errorModelFor = (seq) => stateFor(seq).error;
  /** @param {LinearSeqLike} seq */
  // A stored '' (no error) is no model: it is returned as is.
  const errorFor = (seq) => {
    const error = errorModelFor(seq);
    return error && typeof error === 'object' ? errorText(error.summaryParts) : error;
  };
  /** @param {LinearSeqLike} seq */
  const warningFor = (seq) => {
    const selectorState = stateFor(seq);
    const currentValue = currentSelectorValue(seq);
    if (selectorState.status !== 'ready' || !currentValue) return '';
    const available = buildRecordOptions(selectorState.records).some((option) => option.value === currentValue);
    return available ? '' : t('Selected record was not found in the current file. Choose another record.');
  };

  return {
    selectorStateByUid,
    refresh,
    purgeInactiveState,
    optionsFor,
    statusFor,
    recordsFor,
    isDisabled,
    errorFor,
    errorModelFor,
    warningFor
  };
};
