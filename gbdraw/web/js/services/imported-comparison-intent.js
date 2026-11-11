// @ts-check
/** @import { LinearSeq } from './drawing/inputs.js' */
/** @import { LinearComparisonEdge, LinearComparisonLayoutRow, LinearComparisonResolution } from './linear-comparisons.js' */
/** @import { LinearDrawnUnit, LinearUnitPlan, LinearUnitPlanPayload, LinearUnitRecord, LinearUnitSearch, LinearCardPlacesOf } from './linear-units.js' */
/** @import { CanonicalRecordWindow } from '../contracts/request.generated.js' */
import { liftLegacyRenderRequest } from './legacy-render-request.js';
import { COMPARISON_THRESHOLD_KEYS } from '../mode-profiles.js';
import {
  LINEAR_COMPARISON_MODES, LINEAR_COMPARISON_SOURCES, createLinearComparisonEdge, resolveLinearComparisonPlan
} from './linear-comparisons.js';
import {
  drawsSeveralUnits, planLinearUnitSearch, readLinearUnitPlan, savedPairsUnitSearch, unitRecordInputs, unitRecordsOf
} from './linear-units.js';
import { readStoredCollinearityParams } from './losat-normalization.js';
import { CANONICAL_REQUEST_DESCRIPTION } from '../contracts/request-shape.generated.js';
import { shapeViolations } from '../utils/contract-shape.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { msg } from '../utils/i18n.js';
import { toRaw } from '../utils/vue-raw.js';

/**
 * Whether and how a drawing's saved comparison can be edited, as the Imported comparison state of the drawing holds
 * it. `source` is what Inherit copies (`savedComparisonSource` of the request the comparison was classified on); a
 * comparison that is not preserved has none.
 * @typedef {{
 *   disposition: string,
 *   action: string | null,
 *   message: string,
 *   hasCommittedComparison: boolean,
 *   source: SavedComparisonSource | null
 * }} ImportedComparisonIntent
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
/**
 * One saved record: its key, its row, and its windows (a request of schema 9 or older wrote none: the record draws
 * one unit, as Python reads it).
 * @typedef {{ recordKey?: unknown, presentation?: { gridRow?: unknown } | null, windows?: unknown }} SavedRecord
 */
/**
 * The request records of one saved endpoint (`queryRecordIndex`, `subjectRecordIndex`, which count drawn units),
 * by index into the request's records.
 * @typedef {{ query: number, subject: number }} RecordPair
 */
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
/** @typedef {{ recordKey: string, windows: unknown[] }} SavedUnitInput A record as Python plans its drawn units. */
/**
 * What Inherit copies into a Generate's request, read from the request a drawing's comparison was classified on: its
 * comparisons, the records (keys and windows) whose drawn units their indexes count, the score thresholds it sets
 * (a field it leaves out is left out), and the resources its comparisons read, by ID.
 * @typedef {{
 *   comparisons: SavedComparison[],
 *   records: SavedUnitInput[],
 *   thresholds: { [field: string]: unknown },
 *   resources: { [resourceId: string]: unknown }
 * }} SavedComparisonSource
 */
/**
 * What Generate's Inherit searches for a saved LOSATP pipeline (`inheritedPipelineSearch`): the pipeline, the searches
 * of its unit pairs, and the score thresholds the saved request sets (`SavedComparisonSource.thresholds`).
 * `units` are Python's plans of the saved and the drawn units it carried the pairs by, which the request's Inherit
 * carries the saved comparison's other indexes by too (`inheritSavedComparison`).
 * @typedef {{ pipeline: SavedPipeline, unitSearch: LinearUnitSearch, thresholds: SavedComparisonSource['thresholds'],
 *   units: InheritUnitPlans | null }} InheritedPipelineSearch
 */
/** Runs the Worker helper `planLinearUnits`. @typedef {(payload: LinearUnitPlanPayload) => Promise<{ result?: unknown }>} LinearUnitPlanner */
/**
 * The drawn units Python plans (`planLinearUnits`) for the saved records and for the records a request draws: Inherit
 * carries each saved unit index to the drawn unit of the same unit key and window (`inheritUnitPlans`).
 * @typedef {{ saved: LinearUnitPlan, candidate: LinearUnitPlan }} InheritUnitPlans
 */
/**
 * A History checkpoint of a commit that carries the saved comparison: the caller's checkpoint and the drawn units of
 * what Inherit copied at that point (`withComparisonSourceCheckpoint`, which keeps the source itself apart).
 * @typedef {{ caller: unknown, records: SavedUnitInput[] | null }} ComparisonSourceCheckpoint
 */
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
  hasCommittedComparison: false,
  source: null
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
 * The records a saved endpoint compares: its indexes count drawn units (one per window, D-01), each read as the
 * request record it draws (`unitRecordsOf`). Without the request's units (`unitRecords` null: a record is drawn in
 * several windows and Python's units are not known) an index is only checked to be a unit index, and the records are
 * unknown (`pair` null).
 * @param {SavedEndpoint | null | undefined} comparison
 * @param {(LinearUnitRecord | null)[] | null} unitRecords
 * @returns {{ issue: string } | { pair: RecordPair | null }}
 */
const endpointRecords = (comparison, unitRecords) => {
  const query = Number(comparison?.queryRecordIndex);
  const subject = Number(comparison?.subjectRecordIndex);
  const issue = { issue: msg('The saved comparison does not identify two available records.') };
  if (!Number.isInteger(query) || !Number.isInteger(subject) || query < 0 || subject < 0 || query === subject) return issue;
  if (unitRecords === null) return { pair: null };
  const queryRecord = unitRecords[query];
  const subjectRecord = unitRecords[subject];
  if (!queryRecord || !subjectRecord || queryRecord.recordIndex === subjectRecord.recordIndex) return issue;
  return { pair: { query: queryRecord.recordIndex, subject: subjectRecord.recordIndex } };
};

/**
 * @param {SavedComparison} comparison
 * @param {(LinearUnitRecord | null)[] | null} unitRecords The request record each drawn unit draws (`unitRecordsOf`).
 * @param {number} requestSchema
 */
const pipelineIssue = (comparison, unitRecords, requestSchema) => {
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
    const read = endpointRecords(pair, unitRecords);
    if ('issue' in read) return read.issue;
  }
  return '';
};

/**
 * @param {SavedRecord[]} records
 * @param {readonly RecordPair[]} pairs
 */
const rowsAreEditable = (records, pairs) => {
  const rows = records.map((record, index) => (
    positiveRow(record?.presentation?.gridRow, index + 1)
  ));
  const orderedRows = [...new Set(rows)].sort((left, right) => left - right);
  const rowPosition = new Map(orderedRows.map((row, index) => [row, index]));
  return pairs.every(({ query, subject }) => (
    Math.abs((rowPosition.get(rows[query]) ?? NaN) - (rowPosition.get(rows[subject]) ?? NaN)) === 1
  ));
};

/**
 * @param {SavedRecord[]} records
 * @param {readonly RecordPair[]} pairs
 */
const endpointKeys = (records, pairs) => pairs.map(({ query, subject }) => (
  endpointKey(recordKey(records[query], query), recordKey(records[subject], subject))
));

/**
 * The request's records as Python plans their drawn units (`planLinearUnits`): each record's key and windows.
 * @param {readonly SavedRecord[]} records
 * @returns {SavedUnitInput[]}
 */
const unitInputsOf = (records) => records.map((record, index) => ({
  recordKey: recordKey(record, index), windows: Array.isArray(record?.windows) ? record.windows : []
}));

/**
 * The LOSATP mode whose results a request stores beside a pipeline marker of mode `none`, as the CLI and the Web
 * write them after a render: Collinear for a `collinearityResult`, Similarity groups for an `orthogroupResult`,
 * Pairwise for protein tables alone, and '' for none. The one reader of that mode: the classifier, Load's draft of
 * a CLI Session, which has no saved controls (`projectCanonicalSessionRequest` in services/session-request.js), and
 * Generate's reuse of resolved protein results (`canReuseResolvedProteinArtifacts` in app/run-analysis.js).
 * @param {readonly { kind?: unknown }[]} comparisons The request's comparisons.
 * @returns {'collinear' | 'orthogroup' | 'pairwise' | ''}
 */
export const storedProteinResultMode = (comparisons) => {
  const stores = (/** @type {string} */ kind) => comparisons.some((comparison) => comparison?.kind === kind);
  if (stores('collinearityResult')) return 'collinear';
  if (stores('orthogroupResult')) return 'orthogroup';
  return stores('precomputedProteinComparison') ? 'pairwise' : '';
};

/** @param {{ comparisons: SavedComparison[], pipeline: SavedPipeline }} saved Its rows are editable. */
const generatedPipelineIsProjectable = ({ comparisons, pipeline }) => {
  const endpoints = comparisons.filter((comparison) => ENDPOINT_KINDS.has(comparison.kind));
  if (pipeline.mode === 'pairwise' && pipeline.pairs.length === 0 && endpoints.length === 0) {
    return false;
  }
  // A disabled pipeline can carry CLI settings without any comparison work.
  if (pipeline.mode === 'none' && !storedProteinResultMode(comparisons)
    && (comparisons.length !== 1 || pipeline.pairs.length !== 0)) return false;
  const collinearity = readStoredCollinearityParams(pipeline.settings?.collinearityParams);
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
 * schema-2 re-save of the same Session does. Its comparison indexes count
 * drawn units (OV-2995), read as records through `unitRecordsOf`.
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
    hasCommittedComparison,
    source: null
  });
  // Only a preserved comparison can be inherited: it holds what Inherit copies.
  /** @param {string} message */
  const preserved = (message) => ({
    disposition: IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY,
    action: null,
    message,
    hasCommittedComparison,
    source: savedComparisonSource({ renderRequest, resources })
  });
  /** @param {string} message */
  const decision = (message) => ({
    disposition: IMPORTED_COMPARISON_DISPOSITIONS.DECISION_REQUIRED,
    action: null,
    message,
    hasCommittedComparison,
    source: null
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

  const unitInputs = unitInputsOf(records);
  // Python's units of the request are not known here: they are the records unless a record draws several windows.
  const unitRecords = unitRecordsOf(unitInputs, null);
  /** @param {SavedEndpoint} endpoint An endpoint the loop below accepted. */
  const recordPairOf = (endpoint) => {
    const read = endpointRecords(endpoint, unitRecords);
    return 'pair' in read ? read.pair : null;
  };

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
      const issue = pipelineIssue(comparison, unitRecords, Number(renderRequest?.schema));
      if (issue) return decision(issue);
      if (pipeline) return decision(msg('The saved comparison contains more than one protein pipeline.'));
      // `pipelineIssue` returned no issue, so the settings and pairs are present.
      pipeline = /** @type {SavedPipeline} */ (comparison);
    } else {
      return decision(`The saved comparison kind '${kind || 'unknown'}' is unsupported.`);
    }
    if (ENDPOINT_KINDS.has(kind)) {
      const read = endpointRecords(comparison, unitRecords);
      if ('issue' in read) return decision(read.issue);
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
  const executablePairs = executableEndpoints.map(recordPairOf);
  if (!rowsAreEditable(records, executablePairs.filter((pair) => pair !== null))) {
    return decision(
      msg('The saved comparison connects records that are not on adjacent rows; the comparison controls support adjacent rows only.')
    );
  }
  // The controls and their projection hold record pairs; the comparisons of a drawing with a record drawn in several
  // windows join drawn units, and run as saved through Inherit.
  if (drawsSeveralUnits(unitInputs) && executableEndpoints.length > 0) {
    return preserved(
      msg('The saved comparison can run, but the current controls cannot edit the comparisons of records drawn in several windows.')
    );
  }
  const duplicateEndpoints = endpointKeys(records,
    /** @type {RecordPair[]} Each unit is a record here: a drawing with several units returned above. */ (endpoints.map(recordPairOf)));
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
      msg('The saved comparison can run, but the current controls cannot edit its stored protein-input format.')
    );
  }
  if (pipeline && !generatedPipelineIsProjectable({ comparisons, pipeline })) {
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
 * Whether the next Generate draws the comparison of the current controls rather than the saved one: an editable
 * comparison, or a saved one the user chose to Replace. Generate then takes the cards from the current controls too,
 * splitting a whole multi-record card into one card per record (OV-3324); a saved comparison it keeps (Inherit) or
 * has not been told what to do with keeps the cards it was saved over, which its record indexes name.
 * @param {ImportedComparisonIntent} intent
 */
export const comparisonFollowsCurrentControls = (intent) => (
  intent.disposition === IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE || intent.action === IMPORTED_COMPARISON_ACTIONS.REPLACE
);

/**
 * @param {{ intent: ImportedComparisonIntent | null | undefined, draftResolution: ComparisonDraftResolution | null | undefined }} state
 */
export const importedComparisonExecution = ({ intent, draftResolution }) => {
  if (intent && comparisonFollowsCurrentControls(intent)) {
    // After Replace, an invalid draft (issues need an edge, so it has an intent) fails in the run
    // with its own diagnostic, as an editable comparison does. An empty draft
    // would clear the saved comparison instead.
    if (intent.action !== IMPORTED_COMPARISON_ACTIONS.REPLACE || draftResolution?.hasComparisonIntent === true) {
      return { ok: true, mode: 'draft' };
    }
    return {
      ok: false,
      message: 'Set one complete, valid comparison before generating.',
      fallbackCode: 'COMPARISON_REPLACEMENT_EMPTY'
    };
  }
  if (intent?.action === IMPORTED_COMPARISON_ACTIONS.INHERIT) {
    return { ok: true, mode: 'inherit' };
  }
  if (intent?.action === IMPORTED_COMPARISON_ACTIONS.CLEAR) {
    return { ok: true, mode: 'clear' };
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
 * What Inherit copies from a saved request (`SavedComparisonSource`): the classifier keeps it with a preserved
 * comparison, and Load replaces it with the same of the request as Load read it (the current schema).
 * @param {{ renderRequest?: SavedRenderRequest, resources?: Readonly<Record<string, unknown>> | null }} saved
 * @returns {SavedComparisonSource}
 */
export const savedComparisonSource = ({ renderRequest, resources }) => {
  const comparisons = Array.isArray(renderRequest?.comparisons) ? renderRequest.comparisons : [];
  const options = renderRequest?.diagramOptions ?? {};
  /** @type {{ [resourceId: string]: unknown }} */
  const read = {};
  collectComparisonResourceIds(comparisons).forEach((resourceId) => {
    if (resources && Object.hasOwn(resources, resourceId)) read[resourceId] = resources[resourceId];
  });
  return {
    comparisons,
    records: unitInputsOf(Array.isArray(renderRequest?.records) ? renderRequest.records : []),
    thresholds: Object.fromEntries(COMPARISON_THRESHOLD_KEYS.filter((field) => Object.hasOwn(options, field))
      .map((field) => [field, options[field]])),
    resources: read
  };
};

/**
 * A commit of a candidate of the committed request (`runCommittedCanonicalCandidate` in app/run-analysis.js: Align, a
 * record rotation, Align at a feature) carries the drawing's saved comparison into the request it commits, which
 * draws its records in the windows and orientations that commit wrote to the File cards. What Inherit copies follows
 * that request, as a Save of it would write it, so the next Generate's Inherit draws the comparison in those windows.
 * A drawing whose comparison is not preserved holds no source and keeps none.
 * @param {ImportedComparisonIntent} intent The drawing's Imported comparison state.
 * @param {{ renderRequest?: SavedRenderRequest, resources?: Readonly<Record<string, unknown>> | null }} committed
 */
export const followCommittedComparison = (intent, committed) => {
  if (intent.source !== null) intent.source = savedComparisonSource(committed);
};

/**
 * The History checkpoint hooks of a commit that `followCommittedComparison` follows: the caller's (the File cards'
 * windows and orientations), with what Inherit copies beside them, so Undo and Redo of the commit give it back with
 * those windows. A drawing whose comparison is no longer preserved (a later Generate replaced it) keeps none.
 * History stringifies a checkpoint to compare and count it against its shared limit, so the checkpoint names only the
 * drawn units of the source (`records`, what such a commit changes); the source itself, with the stored tables'
 * bytes its resources hold, is kept beside the checkpoint by identity (`heldSources`), never in its JSON.
 * @param {ImportedComparisonIntent} intent The drawing's Imported comparison state.
 * @param {{
 *   captureIntentCheckpoint: (context: { phase: 'before' | 'after' }) => unknown,
 *   restoreIntentCheckpoint: (checkpoint: unknown) => unknown
 * }} hooks The caller's hooks, as History takes them (`HistoryStepBase` in services/history.js).
 * @returns {{ captureIntentCheckpoint: (context: { phase: 'before' | 'after' }) => Promise<ComparisonSourceCheckpoint>,
 *   restoreIntentCheckpoint: (checkpoint: unknown) => unknown }}
 */
export const withComparisonSourceCheckpoint = (intent, { captureIntentCheckpoint: capture, restoreIntentCheckpoint: restore }) => {
  /** @type {WeakMap<ComparisonSourceCheckpoint, SavedComparisonSource | null>} */
  const heldSources = new WeakMap();
  return {
    captureIntentCheckpoint: async (context) => {
      const source = intent.source;
      /** @type {ComparisonSourceCheckpoint} */
      const checkpoint = { caller: await capture(context), records: source === null ? null : toRaw(source).records };
      heldSources.set(checkpoint, source);
      return checkpoint;
    },
    restoreIntentCheckpoint: (checkpoint) => {
      const held = /** @type {ComparisonSourceCheckpoint} */ (checkpoint);
      if (intent.disposition === IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY && heldSources.has(held)) {
        intent.source = heldSources.get(held) ?? null;
      }
      return restore(held.caller);
    }
  };
};

/**
 * A comparison Inherit cannot carry into this request asks for the choice in the Comparison panel, with the reason the
 * user can act on when there is one.
 * @param {string} [reason]
 */
const inheritChoiceRequired = (reason) => diagnosticError('COMPARISON_CHOICE_REQUIRED', reason ? { reason } : {},
  { operation: 'generate' });

/**
 * Whether each window is a request window (`CanonicalRecordWindow`, the shape Python's decoder reads).
 * @param {readonly unknown[]} windows
 * @returns {windows is CanonicalRecordWindow[]}
 */
const areRecordWindows = (windows) => windows.every((window) => shapeViolations(window,
  { kind: 'ref', name: 'CanonicalRecordWindow' }, CANONICAL_REQUEST_DESCRIPTION.types, 'window').length === 0);

/**
 * Python's plans of the drawn units of the saved records and of `candidateRecords` (`InheritUnitPlans`), or null when
 * neither draws a record in several windows: each unit is then a record (`drawsSeveralUnits`), found by its key.
 * @param {readonly SavedUnitInput[]} savedRecords
 * @param {readonly SavedUnitInput[]} candidateRecords
 * @param {LinearUnitPlanner} runHelper
 * @returns {Promise<InheritUnitPlans | null>}
 */
const planUnitsOf = async (savedRecords, candidateRecords, runHelper) => {
  if (!drawsSeveralUnits(savedRecords) && !drawsSeveralUnits(candidateRecords)) return null;
  /** @param {readonly SavedUnitInput[]} records */
  const plan = async (records) => {
    /** @type {LinearUnitPlanPayload} */
    const payload = { records: records.map(({ recordKey: key, windows }) => {
      // A window that is no request window (a damaged Session) names no unit: the comparison is not carried.
      if (!areRecordWindows(windows)) throw inheritChoiceRequired('INHERIT_UNITS_CHANGED');
      return { recordKey: key, windows, gridRow: null };
    }), pairs: [] };
    return readLinearUnitPlan((await runHelper(payload)).result, payload);
  };
  return { saved: await plan(savedRecords), candidate: await plan(candidateRecords) };
};

/**
 * The plans `inheritSavedComparison` carries a saved comparison into `candidate` by (Generate's Inherit, Gallery
 * publication's rebuild).
 * @param {{ candidate: { renderRequest?: SavedRenderRequest }, source: SavedComparisonSource | null, runHelper: LinearUnitPlanner }} input
 * @returns {Promise<InheritUnitPlans | null>}
 */
export const inheritUnitPlans = ({ candidate, source: held, runHelper }) => {
  const source = toRaw(held);
  const records = Array.isArray(candidate.renderRequest?.records) ? candidate.renderRequest.records : [];
  return source === null ? Promise.resolve(null) : planUnitsOf(source.records, unitInputsOf(records), runHelper);
};

/** @param {LinearDrawnUnit['window']} saved @param {LinearDrawnUnit['window']} drawn */
const sameWindow = (saved, drawn) => (saved === null || drawn === null ? saved === drawn
  : saved.start === drawn.start && saved.end === drawn.end && saved.reverseComplement === drawn.reverseComplement);

/**
 * The saved index of a drawn unit, carried into a request that draws `candidateRecords`. The saved indexes count drawn
 * units (OV-2995), which Python plans (`units`, `planUnitsOf`). Each saved unit is the candidate's unit of the same
 * unit key and window, wherever the request draws it: Load may order the records differently (OV-3365), and two
 * windows of one record may draw the same span (D-18). A request that draws other units asks for the choice. Without
 * a record drawn in several windows a unit is a record, found by its key.
 * @param {SavedComparisonSource} source
 * @param {readonly SavedUnitInput[]} candidateRecords
 * @param {InheritUnitPlans | null} units Python's plans of the two, when either draws a record in several windows.
 * @returns {(value: unknown) => number}
 */
const savedIndexCarrier = (source, candidateRecords, units) => {
  const committedRecords = source.records;
  if (drawsSeveralUnits(committedRecords) || drawsSeveralUnits(candidateRecords)) {
    if (units === null) throw new Error('Inherit needs the unit plans of records drawn in several windows.');
    const planned = [...new Set(units.candidate.units.map((unit) => unit.recordKey))];
    if (JSON.stringify(planned) !== JSON.stringify(candidateRecords.map(({ recordKey: key }) => key))) {
      throw new Error('The unit plan is not of the records the request draws.');
    }
    const candidateIndex = new Map(units.candidate.units.map((unit, index) => [unit.unitKey, index]));
    const carried = units.saved.units.map((unit) => {
      const index = candidateIndex.get(unit.unitKey);
      return index !== undefined && sameWindow(unit.window, units.candidate.units[index].window) ? index : -1;
    });
    if (carried.length !== candidateIndex.size || carried.includes(-1)) throw inheritChoiceRequired('INHERIT_UNITS_CHANGED');
    return (value) => {
      const index = carried[Number(value)];
      if (!Number.isInteger(Number(value)) || index === undefined) throw inheritChoiceRequired('INHERIT_UNITS_CHANGED');
      return index;
    };
  }
  /** @type {Map<string, number>} */
  const candidateIndexByKey = new Map();
  candidateRecords.forEach(({ recordKey: key }, index) => {
    if (candidateIndexByKey.has(key)) throw inheritChoiceRequired();
    candidateIndexByKey.set(key, index);
  });
  return (value) => {
    const sourceIndex = Number(value);
    if (!Number.isInteger(sourceIndex) || !committedRecords[sourceIndex]) throw inheritChoiceRequired('INHERIT_UNITS_CHANGED');
    const candidateIndex = candidateIndexByKey.get(committedRecords[sourceIndex].recordKey);
    if (candidateIndex === undefined) throw inheritChoiceRequired('INHERIT_UNITS_CHANGED');
    return candidateIndex;
  };
};

/**
 * The saved LOSATP Pairwise pipeline that Generate's Inherit runs in the browser (OV-3003). The browser runs LOSAT
 * only for the searches Generate plans, and Python's render has no LOSAT, so this pipeline is searched by the draft
 * path's planner and runner; any other saved comparison is copied as saved.
 * @param {SavedComparisonSource} source
 * @returns {SavedPipeline | null}
 */
const searchedPipelineOf = (source) => /** @type {SavedPipeline | undefined} */ (source.comparisons.find((comparison) => (
  comparison?.kind === 'generatedProteinComparison' && comparison.mode === 'pairwise'
  && isObject(comparison.settings) && Array.isArray(comparison.pairs)
))) ?? null;

/**
 * Generate's Inherit of a saved LOSATP Pairwise pipeline (`searchedPipelineOf`): the pipeline, and the searches over
 * the units Python plans for the drawn cards, which Generate runs with the pipeline's settings through the planner
 * and runner of the draft path. A pipeline that names unit pairs searches exactly those (`savedPairsUnitSearch`), as
 * a selected plan of their records; one that names none searches what Python's Pairwise searches without pairs, the
 * records of adjacent rows (an adjacent plan, `planLinearUnitSearch`; OV-3037). Null when the drawing's saved
 * comparison holds no such pipeline: Inherit copies it as saved (`inheritSavedComparison`).
 * @param {{
 *   source: SavedComparisonSource | null,
 *   cards: readonly LinearSeq[],
 *   placesOf: LinearCardPlacesOf,
 *   runHelper: (payload: LinearUnitPlanPayload) => Promise<{ result?: unknown }>,
 *   requirePlan: (resolution: LinearComparisonResolution) => LinearComparisonResolution,
 *   plan: { sequences: readonly LinearSeq[], recordsOff: Iterable<unknown>, layout: readonly LinearComparisonLayoutRow[] }
 * }} input The drawn cards, in request order, and every card with the OFF records and the row layout; `requirePlan`
 *   checks the plan of the records as Generate checks the draft's (`requireLinearComparisonPlanSnapshot`), so a plan
 *   the drawing no longer allows (its records moved into one row) stops before any search.
 * @returns {Promise<InheritedPipelineSearch | null>}
 */
export const inheritedPipelineSearch = async ({ source: held, cards, placesOf, runHelper, requirePlan, plan }) => {
  const source = toRaw(held);
  const pipeline = source === null ? null : searchedPipelineOf(source);
  if (source === null || pipeline === null) return null;
  const candidate = unitRecordInputs(cards, placesOf);
  const units = await planUnitsOf(source.records, candidate, runHelper);
  const carry = savedIndexCarrier(source, candidate, units);
  /** @param {LinearComparisonEdge[]} edges The selected record pairs; none plans the records of adjacent rows. */
  const losatPlan = (edges) => requirePlan(resolveLinearComparisonPlan({
    ...plan,
    plan: {
      mode: edges.length > 0 ? LINEAR_COMPARISON_MODES.SELECTED : LINEAR_COMPARISON_MODES.ADJACENT,
      defaultSource: LINEAR_COMPARISON_SOURCES.LOSAT, edges
    },
    losatProgram: 'blastp',
    blastpMode: 'pairwise'
  }));
  const pairs = pipeline.pairs.map((pair) => ({
    queryUnit: carry(pair.queryRecordIndex), subjectUnit: carry(pair.subjectRecordIndex)
  }));
  const unitSearch = pairs.length === 0
    ? await planLinearUnitSearch(cards, losatPlan([]), placesOf, runHelper)
    : await savedPairsUnitSearch(cards, pairs, placesOf, runHelper, (recordPairs) => losatPlan(recordPairs.map(
      ({ queryRecordKey, subjectRecordKey }) => createLinearComparisonEdge({
        queryUid: queryRecordKey, subjectUid: subjectRecordKey, source: LINEAR_COMPARISON_SOURCES.LOSAT, included: true
      }))));
  return { pipeline, unitSearch, thresholds: source.thresholds, units };
};

/**
 * Generate's Inherit: puts the saved comparison of a drawing (`source`, its Imported comparison state's) into the
 * request Generate built without one. A comparison it cannot carry into this request (its drawn windows or records
 * changed) asks for the choice in the Comparison panel (`COMPARISON_CHOICE_REQUIRED`). When Generate `searched` the
 * saved pipeline (`inheritedPipelineSearch`), the request keeps what the search wrote in its place. `units` are
 * Python's plans of the saved and the candidate's drawn units (`inheritUnitPlans`).
 * @template {{ renderRequest?: SavedRenderRequest, resources: SavedResources }} C
 * @param {{ candidate: C, source: SavedComparisonSource | null, searched?: boolean, units?: InheritUnitPlans | null }} inherit
 * @returns {C}
 */
export const inheritSavedComparison = ({ candidate, source: held, searched = false, units = null }) => {
  // The drawing state holds the source behind a Vue proxy; a resource descriptor goes into the request as it is.
  const source = toRaw(held);
  if (!isObject(candidate?.renderRequest)) throw inheritChoiceRequired();
  if (source === null || source.comparisons.length === 0) throw inheritChoiceRequired('INHERIT_SOURCE_MISSING');
  const candidateRecords = Array.isArray(candidate.renderRequest.records)
    ? candidate.renderRequest.records
    : [];
  const remapIndex = savedIndexCarrier(source, unitInputsOf(candidateRecords), units);
  const searchedPipeline = searched ? searchedPipelineOf(source) : null;
  const searchedComparisons = searched && Array.isArray(candidate.renderRequest.comparisons)
    ? candidate.renderRequest.comparisons : [];
  const inheritedComparisons = /** @type {SavedComparison[]} */ (structuredClone(
    source.comparisons.filter((comparison) => comparison !== searchedPipeline)));
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
  candidate.renderRequest.comparisons = [...inheritedComparisons, ...searchedComparisons];
  const candidateOptions = candidate.renderRequest.diagramOptions || {};
  COMPARISON_THRESHOLD_KEYS.forEach((field) => {
    if (Object.hasOwn(source.thresholds, field)) {
      candidateOptions[field] = source.thresholds[field];
    } else {
      delete candidateOptions[field];
    }
  });
  candidate.renderRequest.diagramOptions = candidateOptions;

  collectComparisonResourceIds(inheritedComparisons).forEach((resourceId) => {
    const descriptor = source.resources[resourceId];
    if (!isObject(descriptor)) throw inheritChoiceRequired();
    const existing = candidate.resources[resourceId];
    if (existing && existing !== descriptor) throw inheritChoiceRequired();
    candidate.resources[resourceId] = descriptor;
  });
  return candidate;
};
