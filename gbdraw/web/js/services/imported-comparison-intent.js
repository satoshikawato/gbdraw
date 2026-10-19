// @ts-check
import { liftLegacyRenderRequest } from './legacy-render-request.js';
import { msg } from '../utils/i18n.js';

/**
 * Whether and how a Session's saved comparison can be edited, as the Imported comparison state holds it.
 * @typedef {{ disposition: string, action: string | null, message: string, hasCommittedComparison: boolean }} ImportedComparisonIntent
 */
/** @typedef {{ queryRecordIndex?: unknown, subjectRecordIndex?: unknown }} SavedEndpoint */
/**
 * One saved comparison of the canonical request, as these checks read it.
 * @typedef {{
 *   kind?: string,
 *   resourceId?: unknown,
 *   encoding?: string,
 *   valueKind?: string,
 *   queryRecordIndex?: unknown,
 *   subjectRecordIndex?: unknown,
 *   mode?: string,
 *   settings?: { [field: string]: unknown, collinearityParams?: { kind?: string, parameters?: unknown } | null },
 *   pairs?: SavedEndpoint[]
 * }} SavedComparison
 */
/**
 * A generated protein comparison once `pipelineIssue` has accepted it: its settings and pairs are present.
 * @typedef {SavedComparison & {
 *   settings: NonNullable<SavedComparison['settings']>,
 *   pairs: SavedEndpoint[]
 * }} SavedPipeline
 */
/** @typedef {{ recordKey?: unknown, presentation?: { gridRow?: unknown } | null }} SavedRecord */
/**
 * The part of the canonical request the checks read.
 * @typedef {{
 *   schema?: unknown,
 *   mode?: unknown,
 *   comparisons?: SavedComparison[],
 *   records?: SavedRecord[],
 *   diagramOptions?: { [option: string]: unknown, conservationBlastFiles?: unknown }
 * }} SavedRenderRequest
 */
/** @typedef {{ [resourceId: string]: { kind?: unknown } }} SavedResources */
/** @typedef {{ valid?: boolean, hasComparisonIntent?: boolean, error?: string }} ComparisonDraftResolution The current controls' comparison, validated. */

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

/** @type {Set<unknown>} */
const VALID_ACTIONS = new Set(Object.values(IMPORTED_COMPARISON_ACTIONS));
/** @type {Set<string | undefined>} */
const RESOURCE_BACKED_KINDS = new Set([
  'nucleotideBlast',
  'precomputedProteinComparison',
  'orthogroupResult',
  'collinearityResult'
]);
/** @type {Set<string | undefined>} */
const ENDPOINT_KINDS = new Set([
  'nucleotideBlast',
  'precomputedProteinComparison'
]);
const CURRENT_PIPELINE_SETTING_FIELDS = new Set([
  'collinearityParams',
  'collinearityUnitMode',
  'collinearityAnchorMode',
  'collinearitySearchScope',
  'collinearInferOrthogroups',
  'collinearityColorMode',
  'losatpBin',
  'ncbiBlastpBin',
  'losatpThreads',
  'proteinBlastpMaxHits',
  'proteinBlastpCandidateLimit',
  'orthogroupMembershipMode',
  'orthogroupMemberMaxHits',
  'collinearMaxParalogLinksPerOrthogroup'
]);
const LEGACY_PIPELINE_SETTING_FIELDS = new Set([
  ...CURRENT_PIPELINE_SETTING_FIELDS,
  'alignOrthogroupFeature'
]);
/** @type {Readonly<Record<string, string>>} */
const RESOURCE_KIND_BY_COMPARISON_KIND = Object.freeze({
  nucleotideBlast: 'nucleotide-blast',
  precomputedProteinComparison: 'canonical-tsv',
  orthogroupResult: 'orthogroup-result',
  collinearityResult: 'collinearity-result'
});

/**
 * @param {unknown} value
 * @returns {value is Record<string, unknown>}
 */
const isObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
/**
 * @param {unknown} value
 * @param {number} fallback
 */
const positiveRow = (value, fallback) => {
  const row = Number(value);
  return Number.isInteger(row) && row > 0 ? row : fallback;
};
/**
 * @param {SavedRecord | null | undefined} record
 * @param {number} index
 */
const recordKey = (record, index) => String(record?.recordKey || `record-${index + 1}`);
/**
 * @param {string} queryKey
 * @param {string} subjectKey
 */
const endpointKey = (queryKey, subjectKey) => `${queryKey}->${subjectKey}`;

/** @returns {ImportedComparisonIntent} */
export const createImportedComparisonIntentState = () => ({
  disposition: IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE,
  action: null,
  message: '',
  hasCommittedComparison: false
});

/** @param {{ action?: unknown } | null | undefined} intent */
export const serializeImportedComparisonResolution = (intent) => ({
  // `has` is false for a missing intent, so `intent` is set when it is read.
  action: VALID_ACTIONS.has(intent?.action) ? /** @type {{ action?: unknown }} */ (intent).action : null
});

/**
 * @param {SavedComparison} comparison
 * @param {SavedResources | null | undefined} resources
 */
const comparisonResourceIssue = (comparison, resources) => {
  const resourceId = String(comparison?.resourceId || '').trim();
  if (!resourceId || !isObject(resources?.[resourceId])) {
    return msg('The saved comparison is missing a required resource.');
  }
  if (resources[resourceId].kind !== RESOURCE_KIND_BY_COMPARISON_KIND[/** @type {string} */ (comparison.kind)]) {
    return msg('The saved comparison resource has an incompatible kind.');
  }
  if (
    comparison.kind === 'precomputedProteinComparison'
    && comparison.encoding !== 'canonicalTsv'
  ) {
    return msg('The saved protein comparison uses an unsupported encoding.');
  }
  if (
    comparison.kind === 'orthogroupResult'
    && comparison.encoding !== 'canonicalJson'
  ) {
    return msg('The saved Similarity Group result uses an unsupported encoding.');
  }
  if (
    comparison.kind === 'collinearityResult'
    && (
      comparison.encoding !== 'canonicalJson'
      || !(/** @type {(string | undefined)[]} */ (['result', 'blocks'])).includes(comparison.valueKind)
    )
  ) {
    return msg('The saved Collinear result uses an unsupported representation.');
  }
  return '';
};

/**
 * @param {SavedEndpoint | null | undefined} comparison
 * @param {number} recordCount
 */
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
    return msg('The saved comparison does not identify two available records.');
  }
  return '';
};

/**
 * @param {SavedComparison} comparison
 * @param {number} recordCount
 * @param {number} requestSchema
 */
const pipelineIssue = (comparison, recordCount, requestSchema) => {
  if (
    !(/** @type {(string | undefined)[]} */ (['none', 'pairwise', 'orthogroup', 'collinear'])).includes(comparison?.mode)
    || !isObject(comparison?.settings)
    || !Array.isArray(comparison?.pairs)
  ) {
    return msg('The saved protein comparison pipeline is incomplete or unsupported.');
  }
  const acceptedFields = requestSchema >= 8
    ? CURRENT_PIPELINE_SETTING_FIELDS
    : LEGACY_PIPELINE_SETTING_FIELDS;
  const settingKeys = Object.keys(comparison.settings);
  if (
    settingKeys.length !== acceptedFields.size - (Object.prototype.hasOwnProperty.call(comparison.settings, 'collinearInferOrthogroups') ? 0 : 1)
    || settingKeys.some((key) => !acceptedFields.has(key))
  ) {
    return msg('The saved protein comparison pipeline is incomplete or unsupported.');
  }
  const collinearity = comparison.settings.collinearityParams;
  if (collinearity !== null) {
    const parameters = collinearity?.parameters;
    const parameterKeys = isObject(parameters) ? Object.keys(parameters) : [];
    if (
      !isObject(collinearity)
      || !(/** @type {unknown[]} */ (['lossless', 'standard'])).includes(collinearity.kind)
      || parameterKeys.length === 0
    ) {
      return msg('The saved Collinear parameters are incomplete or unsupported.');
    }
  }
  for (const pair of comparison.pairs) {
    const issue = endpointIssue(pair, recordCount);
    if (issue) return issue;
  }
  return '';
};

/**
 * @param {SavedRecord[]} records
 * @param {SavedEndpoint[]} comparisons
 */
const rowsAreEditable = (records, comparisons) => {
  const rows = records.map((record, index) => (
    positiveRow(record?.presentation?.gridRow, index + 1)
  ));
  const orderedRows = [...new Set(rows)].sort((left, right) => left - right);
  const rowPosition = new Map(orderedRows.map((row, index) => [row, index]));
  return comparisons.every((comparison) => {
    const query = Number(comparison.queryRecordIndex);
    const subject = Number(comparison.subjectRecordIndex);
    return Math.abs((rowPosition.get(rows[query]) ?? NaN) - (rowPosition.get(rows[subject]) ?? NaN)) === 1;
  });
};

/**
 * @param {SavedRecord[]} records
 * @param {SavedEndpoint[]} comparisons
 */
const endpointKeys = (records, comparisons) => comparisons.map((comparison) => (
  endpointKey(
    recordKey(records[Number(comparison.queryRecordIndex)], Number(comparison.queryRecordIndex)),
    recordKey(records[Number(comparison.subjectRecordIndex)], Number(comparison.subjectRecordIndex))
  )
));

/** @param {{ records: SavedRecord[], comparisons: SavedComparison[], pipeline: SavedPipeline }} saved */
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
    // A disabled pipeline can carry CLI settings without any comparison work.
    if (!inferredMode && (comparisons.length !== 1 || pipeline.pairs.length !== 0)) return false;
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

/**
 * A stored request is classified as its reader reads it: a schema 1 or 2
 * request through the lift (PR-D, OV-796), so its comparisons classify as the
 * schema-2 re-save of the same Session does.
 * @param {{ renderRequest?: SavedRenderRequest, resources?: SavedResources }} [input]
 * @returns {ImportedComparisonIntent}
 */
export const classifyImportedComparisonIntent = ({
  renderRequest: stored,
  resources = {}
} = {}) => {
  // The cast: a lifted request is a schema-5 request, whose parts the checks
  // below read as they read a stored one.
  const renderRequest = /** @type {SavedRenderRequest | undefined} */ (liftLegacyRenderRequest(stored).value);
  const mode = String(renderRequest?.mode || '');
  const comparisons = Array.isArray(renderRequest?.comparisons)
    ? renderRequest.comparisons
    : [];
  const conservation = Array.isArray(renderRequest?.diagramOptions?.conservationBlastFiles)
    ? renderRequest.diagramOptions.conservationBlastFiles
    : [];
  const hasCommittedComparison = comparisons.length > 0 || conservation.length > 0;
  /** @param {string} [message] */
  const editable = (message = '') => ({
    disposition: IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE,
    action: null,
    message,
    hasCommittedComparison
  });
  /** @param {string} message */
  const preserved = (message) => ({
    disposition: IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY,
    action: null,
    message,
    hasCommittedComparison
  });
  /** @param {string} message */
  const decision = (message) => ({
    disposition: IMPORTED_COMPARISON_DISPOSITIONS.DECISION_REQUIRED,
    action: null,
    message,
    hasCommittedComparison
  });

  if (!hasCommittedComparison) return editable();
  if (!['circular', 'linear'].includes(mode)) {
    return decision(msg('The saved comparison has no supported diagram mode.'));
  }

  if (mode === 'circular') {
    if (comparisons.length > 0) {
      return decision(msg('A Circular session contains unsupported Linear comparison data.'));
    }
    for (const reference of conservation) {
      const resourceId = String(reference?.resourceId || '').trim();
      if (!resourceId || !isObject(resources?.[resourceId])) {
        return decision(msg('The saved Circular comparison is missing a required resource.'));
      }
    }
    return editable(msg('The saved Circular comparison is represented by the current controls.'));
  }

  const records = Array.isArray(renderRequest?.records) ? renderRequest.records : [];
  if (records.length === 0) {
    return decision(msg('The saved comparison has no available records.'));
  }

  /** @type {SavedPipeline | null} */
  let pipeline = null;
  /** @type {Set<string>} */
  const singletonKinds = new Set();
  for (const comparison of comparisons) {
    const kind = String(comparison?.kind || '');
    if (RESOURCE_BACKED_KINDS.has(kind)) {
      const resourceIssue = comparisonResourceIssue(comparison, resources);
      if (resourceIssue) return decision(resourceIssue);
    } else if (kind === 'generatedProteinComparison') {
      const issue = pipelineIssue(comparison, records.length, Number(renderRequest?.schema));
      if (issue) return decision(issue);
      if (pipeline) return decision(msg('The saved comparison contains more than one protein pipeline.'));
      // `pipelineIssue` returned no issue, so the settings and pairs are present.
      pipeline = /** @type {SavedPipeline} */ (comparison);
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
      msg('The saved comparison endpoints are not executable in the supported adjacent-row layout.')
    );
  }
  const duplicateEndpoints = endpointKeys(records, endpoints);
  if (new Set(duplicateEndpoints).size !== duplicateEndpoints.length) {
    return preserved(
      msg('The saved comparison contains repeated record pairs that current controls cannot represent losslessly.')
    );
  }
  const hasDirectProteinInput = comparisons.some((comparison) => (
    (/** @type {(string | undefined)[]} */ (['precomputedProteinComparison', 'orthogroupResult', 'collinearityResult']))
      .includes(comparison.kind)
  ));
  if (hasDirectProteinInput && !pipeline) {
    return preserved(
      msg('The saved comparison is executable, but current controls cannot edit this canonical protein input.')
    );
  }
  if (pipeline && !generatedPipelineIsProjectable({ records, comparisons, pipeline })) {
    return preserved(
      msg('The saved protein comparison is executable, but current controls cannot represent it losslessly.')
    );
  }

  return editable(msg('The saved comparison is represented by the current controls.'));
};

/**
 * @param {ImportedComparisonIntent} target
 * @param {Partial<ImportedComparisonIntent> | null | undefined} classification
 * @param {{ action?: unknown } | null} [storedResolution]
 * @returns {ImportedComparisonIntent}
 */
export const restoreImportedComparisonIntent = (
  target,
  classification,
  storedResolution = /** @type {{ action?: unknown } | null} */ (null)
) => {
  Object.assign(target, createImportedComparisonIntentState(), classification || {});
  const storedAction = String(storedResolution?.action || '').trim().toUpperCase();
  /** @type {Set<string>} */
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

/**
 * @param {{ intent: ImportedComparisonIntent, action: unknown, draftResolution: ComparisonDraftResolution | null | undefined }} choice
 * @returns {{ ok: false, message: string } | { ok: true, action: string }}
 */
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

/**
 * @param {{ intent: ImportedComparisonIntent | null | undefined, draftResolution: ComparisonDraftResolution | null | undefined }} state
 */
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
    // An invalid draft (issues need an edge, so it has an intent) fails in the run
    // with its own diagnostic, as an editable comparison does. An empty draft
    // would clear the saved comparison instead.
    if (draftResolution?.hasComparisonIntent === true) {
      return { ok: true, mode: 'draft' };
    }
    return {
      ok: false,
      message: 'Set one complete, valid comparison before generating.',
      fallbackCode: 'COMPARISON_REPLACEMENT_EMPTY'
    };
  }
  // A reason with its own diagnostic (a missing resource) keeps it; any other
  // reason reports the pending choice instead of the unclassified fallback (OV-303).
  return {
    ok: false,
    message: intent?.message || 'Choose how to resolve the saved comparison before generating.',
    fallbackCode: 'COMPARISON_CHOICE_REQUIRED'
  };
};

/**
 * @param {unknown} value
 * @param {Set<string>} [target]
 * @returns {Set<string>}
 */
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

/**
 * @template {{ renderRequest?: SavedRenderRequest, resources: SavedResources }} C
 * @param {{
 *   candidate: C,
 *   committed: { renderRequest?: SavedRenderRequest, resources?: Readonly<Record<string, unknown>> } | null | undefined
 * }} drafts The committed Session's resources are read by id, and each entry read is checked.
 * @returns {C}
 */
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
  /** @type {Map<string, number>} */
  const candidateIndexByKey = new Map();
  candidateRecords.forEach((record, index) => {
    const key = recordKey(record, index);
    if (candidateIndexByKey.has(key)) {
      throw new Error('The current draft has duplicate record identities. Choose Replace or Clear.');
    }
    candidateIndexByKey.set(key, index);
  });
  /** @param {unknown} value */
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
    const existing = candidate.resources[resourceId];
    if (existing && existing !== descriptor) {
      throw new Error(
        'The preserved comparison resource conflicts with the current draft. Choose Replace or Clear.'
      );
    }
    candidate.resources[resourceId] = descriptor;
  });
  return candidate;
};
