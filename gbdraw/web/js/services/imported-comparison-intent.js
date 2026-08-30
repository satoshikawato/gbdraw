export const IMPORTED_COMPARISON_DISPOSITIONS = Object.freeze({
  EDITABLE: 'EDITABLE',
  PRESERVED_READ_ONLY: 'PRESERVED_READ_ONLY',
  DECISION_REQUIRED: 'DECISION_REQUIRED'
});

export const IMPORTED_COMPARISON_ACTIONS = Object.freeze({
  INHERIT: 'INHERIT',
  REPLACE: 'REPLACE',
  CLEAR: 'CLEAR'
});

const VALID_ACTIONS = new Set(Object.values(IMPORTED_COMPARISON_ACTIONS));
const RESOURCE_BACKED_KINDS = new Set([
  'nucleotideBlast',
  'precomputedProteinComparison',
  'orthogroupResult',
  'collinearityResult'
]);
const ENDPOINT_KINDS = new Set([
  'nucleotideBlast',
  'precomputedProteinComparison'
]);
const PIPELINE_SETTING_FIELDS = new Set([
  'collinearityParams',
  'collinearityUnitMode',
  'collinearityAnchorMode',
  'collinearitySearchScope',
  'collinearityColorMode',
  'losatpBin',
  'ncbiBlastpBin',
  'losatpThreads',
  'proteinBlastpMaxHits',
  'proteinBlastpCandidateLimit',
  'orthogroupMembershipMode',
  'orthogroupMemberMaxHits',
  'collinearMaxParalogLinksPerOrthogroup',
  'alignOrthogroupFeature'
]);
const RESOURCE_KIND_BY_COMPARISON_KIND = Object.freeze({
  nucleotideBlast: 'nucleotide-blast',
  precomputedProteinComparison: 'canonical-tsv',
  orthogroupResult: 'orthogroup-result',
  collinearityResult: 'collinearity-result'
});

const isObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
const positiveRow = (value, fallback) => {
  const row = Number(value);
  return Number.isInteger(row) && row > 0 ? row : fallback;
};
const recordKey = (record, index) => String(record?.recordKey || `record-${index + 1}`);
const endpointKey = (queryKey, subjectKey) => `${queryKey}->${subjectKey}`;

export const createImportedComparisonIntentState = () => ({
  disposition: IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE,
  action: null,
  message: '',
  hasCommittedComparison: false
});

export const serializeImportedComparisonResolution = (intent) => ({
  action: VALID_ACTIONS.has(intent?.action) ? intent.action : null
});

const comparisonResourceIssue = (comparison, resources) => {
  const resourceId = String(comparison?.resourceId || '').trim();
  if (!resourceId || !isObject(resources?.[resourceId])) {
    return 'The saved comparison is missing a required resource.';
  }
  if (resources[resourceId].kind !== RESOURCE_KIND_BY_COMPARISON_KIND[comparison.kind]) {
    return 'The saved comparison resource has an incompatible kind.';
  }
  if (
    comparison.kind === 'precomputedProteinComparison'
    && comparison.encoding !== 'canonicalTsv'
  ) {
    return 'The saved protein comparison uses an unsupported encoding.';
  }
  if (
    comparison.kind === 'orthogroupResult'
    && comparison.encoding !== 'canonicalJson'
  ) {
    return 'The saved Similarity Group result uses an unsupported encoding.';
  }
  if (
    comparison.kind === 'collinearityResult'
    && (
      comparison.encoding !== 'canonicalJson'
      || !['result', 'blocks'].includes(comparison.valueKind)
    )
  ) {
    return 'The saved Collinear result uses an unsupported representation.';
  }
  return '';
};

const endpointIssue = (comparison, recordCount) => {
  const query = Number(comparison?.queryRecordIndex);
  const subject = Number(comparison?.subjectRecordIndex);
  if (
    !Number.isInteger(query)
    || !Number.isInteger(subject)
    || query < 0
    || subject < 0
    || query >= recordCount
    || subject >= recordCount
    || query === subject
  ) {
    return 'The saved comparison does not identify two available records.';
  }
  return '';
};

const pipelineIssue = (comparison, recordCount) => {
  if (
    !['none', 'pairwise', 'orthogroup', 'collinear'].includes(comparison?.mode)
    || !isObject(comparison?.settings)
    || !Array.isArray(comparison?.pairs)
  ) {
    return 'The saved protein comparison pipeline is incomplete or unsupported.';
  }
  const settingKeys = Object.keys(comparison.settings);
  if (
    settingKeys.length !== PIPELINE_SETTING_FIELDS.size
    || settingKeys.some((key) => !PIPELINE_SETTING_FIELDS.has(key))
  ) {
    return 'The saved protein comparison pipeline is incomplete or unsupported.';
  }
  const collinearity = comparison.settings.collinearityParams;
  if (collinearity !== null) {
    const parameters = collinearity?.parameters;
    const parameterKeys = isObject(parameters) ? Object.keys(parameters) : [];
    if (
      !isObject(collinearity)
      || !['lossless', 'standard'].includes(collinearity.kind)
      || parameterKeys.length === 0
    ) {
      return 'The saved Collinear parameters are incomplete or unsupported.';
    }
  }
  for (const pair of comparison.pairs) {
    const issue = endpointIssue(pair, recordCount);
    if (issue) return issue;
  }
  return '';
};

const rowsAreEditable = (records, comparisons) => {
  const rows = records.map((record, index) => (
    positiveRow(record?.presentation?.gridRow, index + 1)
  ));
  const orderedRows = [...new Set(rows)].sort((left, right) => left - right);
  const rowPosition = new Map(orderedRows.map((row, index) => [row, index]));
  return comparisons.every((comparison) => {
    const query = Number(comparison.queryRecordIndex);
    const subject = Number(comparison.subjectRecordIndex);
    return Math.abs(rowPosition.get(rows[query]) - rowPosition.get(rows[subject])) === 1;
  });
};

const endpointKeys = (records, comparisons) => comparisons.map((comparison) => (
  endpointKey(
    recordKey(records[Number(comparison.queryRecordIndex)], Number(comparison.queryRecordIndex)),
    recordKey(records[Number(comparison.subjectRecordIndex)], Number(comparison.subjectRecordIndex))
  )
));

const generatedPipelineIsProjectable = ({ records, comparisons, pipeline }) => {
  const endpoints = comparisons.filter((comparison) => ENDPOINT_KINDS.has(comparison.kind));
  if (!rowsAreEditable(records, [...endpoints, ...pipeline.pairs])) return false;
  if (pipeline.mode === 'pairwise' && pipeline.pairs.length === 0 && endpoints.length === 0) {
    return false;
  }
  if (pipeline.mode === 'none') {
    const inferredMode = comparisons.some((comparison) => comparison.kind === 'collinearityResult')
      ? 'collinear'
      : comparisons.some((comparison) => comparison.kind === 'orthogroupResult')
        ? 'orthogroup'
        : endpoints.some((comparison) => comparison.kind === 'precomputedProteinComparison')
          ? 'pairwise'
          : '';
    if (!inferredMode) return false;
  }
  const collinearity = pipeline.settings?.collinearityParams;
  if (!isObject(collinearity) || collinearity.kind !== 'lossless') return false;
  const parameters = collinearity.parameters;
  if (
    !isObject(parameters)
    || !['minAnchors', 'maxUnitGap', 'maxDiagonalDrift', 'maxConflicts', 'mergeOrientation']
      .every((key) => Object.prototype.hasOwnProperty.call(parameters, key))
  ) return false;
  return true;
};

export const classifyImportedComparisonIntent = ({
  renderRequest,
  resources = {}
} = {}) => {
  const mode = String(renderRequest?.mode || '');
  const comparisons = Array.isArray(renderRequest?.comparisons)
    ? renderRequest.comparisons
    : [];
  const conservation = Array.isArray(renderRequest?.diagramOptions?.conservationBlastFiles)
    ? renderRequest.diagramOptions.conservationBlastFiles
    : [];
  const hasCommittedComparison = comparisons.length > 0 || conservation.length > 0;
  const editable = (message = '') => ({
    disposition: IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE,
    action: null,
    message,
    hasCommittedComparison
  });
  const preserved = (message) => ({
    disposition: IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY,
    action: null,
    message,
    hasCommittedComparison
  });
  const decision = (message) => ({
    disposition: IMPORTED_COMPARISON_DISPOSITIONS.DECISION_REQUIRED,
    action: null,
    message,
    hasCommittedComparison
  });

  if (!hasCommittedComparison) return editable();
  if (!['circular', 'linear'].includes(mode)) {
    return decision('The saved comparison has no supported diagram mode.');
  }

  if (mode === 'circular') {
    if (comparisons.length > 0) {
      return decision('A Circular session contains unsupported Linear comparison data.');
    }
    for (const reference of conservation) {
      const resourceId = String(reference?.resourceId || '').trim();
      if (!resourceId || !isObject(resources?.[resourceId])) {
        return decision('The saved Circular comparison is missing a required resource.');
      }
    }
    return editable('The saved Circular comparison is represented by the current controls.');
  }

  const records = Array.isArray(renderRequest?.records) ? renderRequest.records : [];
  if (records.length === 0) {
    return decision('The saved comparison has no available records.');
  }

  let pipeline = null;
  const singletonKinds = new Set();
  for (const comparison of comparisons) {
    const kind = String(comparison?.kind || '');
    if (RESOURCE_BACKED_KINDS.has(kind)) {
      const resourceIssue = comparisonResourceIssue(comparison, resources);
      if (resourceIssue) return decision(resourceIssue);
    } else if (kind === 'generatedProteinComparison') {
      const issue = pipelineIssue(comparison, records.length);
      if (issue) return decision(issue);
      if (pipeline) return decision('The saved comparison contains more than one protein pipeline.');
      pipeline = comparison;
    } else {
      return decision(`The saved comparison kind '${kind || 'unknown'}' is unsupported.`);
    }
    if (ENDPOINT_KINDS.has(kind)) {
      const issue = endpointIssue(comparison, records.length);
      if (issue) return decision(issue);
    }
    if (['orthogroupResult', 'collinearityResult'].includes(kind)) {
      if (singletonKinds.has(kind)) {
        return decision(`The saved comparison contains more than one ${kind}.`);
      }
      singletonKinds.add(kind);
    }
  }

  const endpoints = comparisons.filter((comparison) => ENDPOINT_KINDS.has(comparison.kind));
  const executableEndpoints = pipeline
    ? [...endpoints, ...pipeline.pairs]
    : endpoints;
  if (!rowsAreEditable(records, executableEndpoints)) {
    return decision(
      'The saved comparison endpoints are not executable in the supported adjacent-row layout.'
    );
  }
  const duplicateEndpoints = endpointKeys(records, endpoints);
  if (new Set(duplicateEndpoints).size !== duplicateEndpoints.length) {
    return preserved(
      'The saved comparison contains repeated record pairs that current controls cannot represent losslessly.'
    );
  }
  const hasDirectProteinInput = comparisons.some((comparison) => (
    ['precomputedProteinComparison', 'orthogroupResult', 'collinearityResult']
      .includes(comparison.kind)
  ));
  if (hasDirectProteinInput && !pipeline) {
    return preserved(
      'The saved comparison is executable, but current controls cannot edit this canonical protein input.'
    );
  }
  if (pipeline && !generatedPipelineIsProjectable({ records, comparisons, pipeline })) {
    return preserved(
      'The saved protein comparison is executable, but current controls cannot represent it losslessly.'
    );
  }

  return editable('The saved comparison is represented by the current controls.');
};

export const restoreImportedComparisonIntent = (
  target,
  classification,
  storedResolution = null
) => {
  Object.assign(target, createImportedComparisonIntentState(), classification || {});
  const storedAction = String(storedResolution?.action || '').trim().toUpperCase();
  const allowed = target.disposition === IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY
    ? new Set(Object.values(IMPORTED_COMPARISON_ACTIONS))
    : target.disposition === IMPORTED_COMPARISON_DISPOSITIONS.DECISION_REQUIRED
      ? new Set([
          IMPORTED_COMPARISON_ACTIONS.REPLACE,
          IMPORTED_COMPARISON_ACTIONS.CLEAR
        ])
      : new Set();
  target.action = allowed.has(storedAction) ? storedAction : null;
  return target;
};

export const resolveImportedComparisonAction = ({
  intent,
  action,
  draftResolution
}) => {
  const requested = String(action || '').trim().toUpperCase();
  if (!VALID_ACTIONS.has(requested)) {
    return { ok: false, message: 'Choose how to resolve the saved comparison.' };
  }
  if (
    requested === IMPORTED_COMPARISON_ACTIONS.INHERIT
    && (
      intent.disposition !== IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY
      || intent.hasCommittedComparison !== true
    )
  ) {
    return { ok: false, message: 'Only a complete preserved comparison can be inherited.' };
  }
  if (
    requested === IMPORTED_COMPARISON_ACTIONS.REPLACE
    && (draftResolution?.valid !== true || draftResolution?.hasComparisonIntent !== true)
  ) {
    return {
      ok: false,
      message: draftResolution?.error
        || 'Set one complete, valid comparison in the current controls before replacing the saved comparison.'
    };
  }
  if (
    intent.disposition === IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE
    && requested !== IMPORTED_COMPARISON_ACTIONS.CLEAR
  ) {
    return { ok: false, message: 'The saved comparison is already editable.' };
  }
  intent.action = requested;
  return { ok: true, action: requested };
};

export const importedComparisonExecution = ({ intent, draftResolution }) => {
  if (intent?.disposition === IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE) {
    return { ok: true, mode: 'draft' };
  }
  if (intent?.action === IMPORTED_COMPARISON_ACTIONS.INHERIT) {
    return { ok: true, mode: 'inherit' };
  }
  if (intent?.action === IMPORTED_COMPARISON_ACTIONS.CLEAR) {
    return { ok: true, mode: 'clear' };
  }
  if (intent?.action === IMPORTED_COMPARISON_ACTIONS.REPLACE) {
    if (draftResolution?.valid === true && draftResolution?.hasComparisonIntent === true) {
      return { ok: true, mode: 'draft' };
    }
    return {
      ok: false,
      message: draftResolution?.error
        || 'Set one complete, valid comparison before generating.'
    };
  }
  return {
    ok: false,
    message: intent?.message || 'Choose how to resolve the saved comparison before generating.'
  };
};

const collectComparisonResourceIds = (value, target = new Set()) => {
  if (Array.isArray(value)) {
    value.forEach((entry) => collectComparisonResourceIds(entry, target));
    return target;
  }
  if (!isObject(value)) return target;
  if (Object.prototype.hasOwnProperty.call(value, 'resourceId')) {
    const resourceId = String(value.resourceId || '').trim();
    if (resourceId) target.add(resourceId);
  }
  Object.values(value).forEach((entry) => collectComparisonResourceIds(entry, target));
  return target;
};

export const inheritCommittedComparisonIntent = ({ candidate, committed }) => {
  if (!isObject(candidate?.renderRequest) || !isObject(committed?.renderRequest)) {
    throw new Error('The preserved comparison is unavailable. Reload the Session and try again.');
  }
  const comparisons = committed.renderRequest.comparisons;
  if (!Array.isArray(comparisons) || comparisons.length === 0) {
    throw new Error('The preserved comparison is unavailable. Choose Replace or Clear.');
  }
  const committedRecords = Array.isArray(committed.renderRequest.records)
    ? committed.renderRequest.records
    : [];
  const candidateRecords = Array.isArray(candidate.renderRequest.records)
    ? candidate.renderRequest.records
    : [];
  const candidateIndexByKey = new Map();
  candidateRecords.forEach((record, index) => {
    const key = recordKey(record, index);
    if (candidateIndexByKey.has(key)) {
      throw new Error('The current draft has duplicate record identities. Choose Replace or Clear.');
    }
    candidateIndexByKey.set(key, index);
  });
  const remapIndex = (value) => {
    const sourceIndex = Number(value);
    if (!Number.isInteger(sourceIndex) || !committedRecords[sourceIndex]) {
      throw new Error('The preserved comparison references an unavailable record.');
    }
    const key = recordKey(committedRecords[sourceIndex], sourceIndex);
    const candidateIndex = candidateIndexByKey.get(key);
    if (!Number.isInteger(candidateIndex)) {
      throw new Error('The current draft no longer contains every record used by the preserved comparison.');
    }
    return candidateIndex;
  };
  const inheritedComparisons = structuredClone(comparisons);
  inheritedComparisons.forEach((comparison) => {
    if (ENDPOINT_KINDS.has(comparison?.kind)) {
      comparison.queryRecordIndex = remapIndex(comparison.queryRecordIndex);
      comparison.subjectRecordIndex = remapIndex(comparison.subjectRecordIndex);
    }
    if (comparison?.kind === 'generatedProteinComparison') {
      comparison.pairs = (Array.isArray(comparison.pairs) ? comparison.pairs : [])
        .map((pair) => ({
          queryRecordIndex: remapIndex(pair?.queryRecordIndex),
          subjectRecordIndex: remapIndex(pair?.subjectRecordIndex)
        }));
    }
  });
  candidate.renderRequest.comparisons = inheritedComparisons;
  const committedOptions = committed.renderRequest.diagramOptions || {};
  const candidateOptions = candidate.renderRequest.diagramOptions || {};
  ['evalue', 'bitscore', 'identity', 'alignmentLength'].forEach((field) => {
    if (Object.prototype.hasOwnProperty.call(committedOptions, field)) {
      candidateOptions[field] = committedOptions[field];
    } else {
      delete candidateOptions[field];
    }
  });
  candidate.renderRequest.diagramOptions = candidateOptions;

  collectComparisonResourceIds(inheritedComparisons).forEach((resourceId) => {
    const descriptor = committed.resources?.[resourceId];
    if (!isObject(descriptor)) {
      throw new Error('The preserved comparison is missing a required resource.');
    }
    const existing = candidate.resources?.[resourceId];
    if (existing && existing !== descriptor) {
      throw new Error(
        'The preserved comparison resource conflicts with the current draft. Choose Replace or Clear.'
      );
    }
    candidate.resources[resourceId] = descriptor;
  });
  return candidate;
};
