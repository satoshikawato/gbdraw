export const AUTOMATIC_RECORD_OPTION_LABEL = 'Automatic (no explicit selector)';
export const RECORDS_LOADING_LABEL = 'Loading records...';
export const RECORD_FILE_REQUIRED_LABEL = 'Upload a sequence file first';
export const RECORD_READ_ERROR_LABEL = 'Records could not be loaded';

const currentSelectorValue = (seq) => String(seq?.region_record_id ?? '').trim();

export const formatRecordLength = (value) => {
  const numeric = Number(value);
  if (!Number.isInteger(numeric) || numeric <= 0) return 'length unavailable';
  return `${numeric.toLocaleString('en-US')} bp`;
};

export const buildRecordOptions = (records, currentValue = '') => {
  const normalizedRecords = Array.isArray(records) ? records : [];
  const idCounts = new Map();
  normalizedRecords.forEach((record) => {
    const recordId = String(record?.recordId ?? '').trim();
    idCounts.set(recordId, (idCounts.get(recordId) || 0) + 1);
  });

  const options = normalizedRecords.map((record, index) => {
    const selector = String(record?.selector ?? `#${index + 1}`).trim() || `#${index + 1}`;
    const recordId = String(record?.recordId ?? '').trim() || `Record_${index + 1}`;
    const duplicate = (idCounts.get(recordId) || 0) > 1;
    return {
      value: duplicate ? selector : recordId,
      label: `${recordId} (${formatRecordLength(record?.recordLength)})${duplicate ? ` [${selector}]` : ''}`,
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
  error: ''
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
  const sourceFileFor = (seq, inputType = state.lInputType.value) =>
    inputType === 'gff' ? seq?.fasta : seq?.gb;

  const stateFor = (seq) => {
    const uid = uidFor(seq);
    return (uid && selectorStateByUid[uid]) || emptySelectorState();
  };

  const replaceState = (uid, nextState) => {
    selectorStateByUid[uid] = {
      status: nextState.status,
      records: Array.isArray(nextState.records) ? nextState.records : [],
      error: String(nextState.error || '')
    };
  };

  const purgeInactiveState = () => {
    const activeUids = new Set(state.linearSeqs.map(uidFor).filter(Boolean));
    Object.keys(selectorStateByUid).forEach((uid) => {
      if (!activeUids.has(uid)) delete selectorStateByUid[uid];
    });
  };

  const isCurrentRequest = ({ generation, uid, file, inputType }) => {
    if (generation !== refreshGeneration) return false;
    if (state.mode.value !== 'linear' || state.lInputType.value !== inputType) return false;
    const currentSeq = state.linearSeqs.find((seq) => uidFor(seq) === uid);
    return Boolean(currentSeq) && sourceFileFor(currentSeq, inputType) === file;
  };

  const refresh = async () => {
    const generation = ++refreshGeneration;
    purgeInactiveState();

    if (state.mode.value !== 'linear') {
      Object.keys(selectorStateByUid).forEach((uid) => delete selectorStateByUid[uid]);
      return;
    }

    const inputType = state.lInputType.value;
    const format = inputType === 'gff' ? 'fasta' : 'genbank';
    const targets = [];
    for (const seq of state.linearSeqs) {
      const uid = uidFor(seq);
      if (!uid) continue;
      const file = sourceFileFor(seq, inputType);
      if (!file) {
        replaceState(uid, emptySelectorState());
        continue;
      }
      replaceState(uid, { status: 'loading', records: [], error: '' });
      targets.push({ uid, file });
    }
    if (!state.pyodideReady.value) return;

    for (const { uid, file } of targets) {
      try {
        const records = await recordReader({
          file,
          format,
          temporaryPath: `/record-selector-${sanitizePathSegment(uid)}-${generation}.${format === 'fasta' ? 'fasta' : 'gb'}`
        });
        if (!isCurrentRequest({ generation, uid, file, inputType })) return;
        replaceState(uid, { status: 'ready', records, error: '' });
      } catch (error) {
        if (!isCurrentRequest({ generation, uid, file, inputType })) return;
        logger.warn?.(`Failed to read records for ${uid}:`, error);
        replaceState(uid, {
          status: 'error',
          records: [],
          error: 'Could not read records from this file.'
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
    isDisabled,
    errorFor,
    warningFor
  };
};
