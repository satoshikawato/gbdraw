import { buildDisambiguatedRecordEntries, formatRecordLength } from './record-options.js';

export const AUTOMATIC_RECORD_OPTION_LABEL = 'Automatic (no explicit selector)';
export const RECORDS_LOADING_LABEL = 'Loading records...';
export const RECORD_FILE_REQUIRED_LABEL = 'Upload a sequence file first';
export const RECORD_READ_ERROR_LABEL = 'Records could not be loaded';

const currentSelectorValue = (seq) => String(seq?.region_record_id ?? '').trim();

export { formatRecordLength };

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
    { value: '', label: AUTOMATIC_RECORD_OPTION_LABEL, synthetic: false },
    ...(selectedIsMissing
      ? [{
          value: selectedValue,
          label: `${selectedValue} (not found in current file)`,
          synthetic: true
        }]
      : []),
    ...options
  ];
};

const emptySelectorState = () => ({
  status: 'idle',
  records: [],
  error: '',
  inputType: '',
  primaryFile: null,
  pairedFile: null
});

const sanitizePathSegment = (value) =>
  String(value || 'sequence').replace(/[^A-Za-z0-9_.-]+/g, '_').slice(0, 96) || 'sequence';

export const createLinearRecordSelector = ({
  state,
  reactive,
  recordReader,
  logger = console
}) => {
  const selectorStateByUid = reactive({});
  let refreshGeneration = 0;

  const uidFor = (seq) => String(seq?.uid ?? '').trim();
  const sourceFilesFor = (seq, inputType = state.lInputType.value) => (
    inputType === 'gff'
      ? { primaryFile: seq?.gff || null, pairedFile: seq?.fasta || null }
      : { primaryFile: seq?.gb || null, pairedFile: null }
  );

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
        status: primaryFile && (inputType !== 'gff' || pairedFile) ? 'loading' : 'idle',
        inputType,
        primaryFile,
        pairedFile
      };
    }
    return stored;
  };

  const replaceState = (uid, nextState) => {
    selectorStateByUid[uid] = {
      status: nextState.status,
      records: Array.isArray(nextState.records) ? nextState.records : [],
      error: String(nextState.error || ''),
      inputType: String(nextState.inputType || ''),
      primaryFile: nextState.primaryFile || null,
      pairedFile: nextState.pairedFile || null
    };
  };

  const purgeInactiveState = () => {
    const activeUids = new Set(state.linearSeqs.map(uidFor).filter(Boolean));
    Object.keys(selectorStateByUid).forEach((uid) => {
      if (!activeUids.has(uid)) delete selectorStateByUid[uid];
    });
  };

  const isCurrentRequest = ({ generation, uid, primaryFile, pairedFile, inputType }) => {
    if (generation !== refreshGeneration) return false;
    if (state.mode.value !== 'linear' || state.lInputType.value !== inputType) return false;
    const currentSeq = state.linearSeqs.find((seq) => uidFor(seq) === uid);
    if (!currentSeq) return false;
    const currentFiles = sourceFilesFor(currentSeq, inputType);
    return currentFiles.primaryFile === primaryFile && currentFiles.pairedFile === pairedFile;
  };

  const refresh = async () => {
    const generation = ++refreshGeneration;
    purgeInactiveState();

    if (state.mode.value !== 'linear') {
      Object.keys(selectorStateByUid).forEach((uid) => delete selectorStateByUid[uid]);
      return;
    }

    const inputType = state.lInputType.value;
    const targets = [];
    for (const seq of state.linearSeqs) {
      const uid = uidFor(seq);
      if (!uid) continue;
      const { primaryFile, pairedFile } = sourceFilesFor(seq, inputType);
      if (!primaryFile || (inputType === 'gff' && !pairedFile)) {
        replaceState(uid, { ...emptySelectorState(), inputType, primaryFile, pairedFile });
        continue;
      }
      replaceState(uid, {
        status: 'loading', records: [], error: '', inputType, primaryFile, pairedFile
      });
      targets.push({ uid, primaryFile, pairedFile });
    }
    if (!state.pyodideReady.value) return;

    for (const { uid, primaryFile, pairedFile } of targets) {
      try {
        const records = await recordReader({
          inputType,
          primaryFile,
          pairedFile,
          temporaryPathPrefix: `/record-selector-${sanitizePathSegment(uid)}-${generation}`
        });
        if (!isCurrentRequest({ generation, uid, primaryFile, pairedFile, inputType })) return;
        replaceState(uid, {
          status: 'ready', records, error: '', inputType, primaryFile, pairedFile
        });
      } catch (error) {
        if (!isCurrentRequest({ generation, uid, primaryFile, pairedFile, inputType })) return;
        logger.warn?.(`Failed to read records for ${uid}:`, error);
        replaceState(uid, {
          status: 'error',
          records: [],
          error: 'Could not read records from this file.',
          inputType,
          primaryFile,
          pairedFile
        });
      }
    }
  };

  const optionsFor = (seq) => {
    const selectorState = stateFor(seq);
    const currentValue = currentSelectorValue(seq);
    if (selectorState.status === 'ready') {
      return buildRecordOptions(selectorState.records, currentValue);
    }

    const label = selectorState.status === 'loading'
      ? RECORDS_LOADING_LABEL
      : selectorState.status === 'error'
        ? RECORD_READ_ERROR_LABEL
        : RECORD_FILE_REQUIRED_LABEL;
    return [{
      value: currentValue,
      label,
      synthetic: Boolean(currentValue)
    }];
  };

  const isDisabled = (seq) => stateFor(seq).status !== 'ready';
  const statusFor = (seq) => stateFor(seq).status;
  const recordsFor = (seq) => stateFor(seq).records.slice();
  const errorFor = (seq) => stateFor(seq).error;
  const warningFor = (seq) => {
    const selectorState = stateFor(seq);
    const currentValue = currentSelectorValue(seq);
    if (selectorState.status !== 'ready' || !currentValue) return '';
    const available = buildRecordOptions(selectorState.records).some((option) => option.value === currentValue);
    return available ? '' : 'Selected record was not found in the current file. Choose another record.';
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
    warningFor
  };
};
