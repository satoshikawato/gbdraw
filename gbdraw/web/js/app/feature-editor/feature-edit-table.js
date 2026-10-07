// @ts-check
import { buildFeatureOverrideTable } from '../run-info.js';
import {
  FEATURE_OVERRIDE_EDIT_FIELDS,
  canonicalFeatureOverrides,
  recordKeyBelongsToRequest,
  requestFeatureOverrides,
  rowBelongsToRequest,
  updateFeatureOverride
} from '../../services/feature-placement.js';
import { diagnosticError, normalizeUserFacingError } from '../../services/error-normalization.js';
import { cloneFileBytesForTransfer } from '../../services/file-content-cache.js';
import { downloadTextFile } from '../../services/text-download.js';

// Export and Load Feature Edits TSV (design Q4 6.4, Owner Q2 = A): the
// per-feature edits of the committed records as a --feature_override_table, so
// the CLI and the Python API read the same file. Python reads a loaded table
// (R4, R9); rows that name no record or feature of the committed records are
// counted and not applied (Owner Q3 = A).
const OPERATION = 'readFeatureOverrideTable';
const CLEARED_EDITS = Object.freeze(Object.fromEntries(FEATURE_OVERRIDE_EDIT_FIELDS.map((field) => [field, null])));
const editsOf = (row) => Object.fromEntries(FEATURE_OVERRIDE_EDIT_FIELDS.map((field) => [field, row[field]]));
const hasEdit = (row) => FEATURE_OVERRIDE_EDIT_FIELDS.some((field) => row?.[field] != null);

// The helper's rows must be request rows of the committed records.
export const admitFeatureOverrideTable = (result, records) => {
  const invalid = () => diagnosticError('RESULT_INVALID', {}, /** @type {{ stage?: string, operation?: string }} */ ({ operation: OPERATION, stage: 'result-admission' }));
  const unmatchedRows = result?.unmatchedRows;
  if (!Array.isArray(result?.rows) || !Array.isArray(unmatchedRows)
    || unmatchedRows.some((row) => !Number.isSafeInteger(row) || row < 2)
    || result.rows.some((row) => !recordKeyBelongsToRequest(row?.recordKey, records))) {
    throw invalid();
  }
  try {
    return { rows: canonicalFeatureOverrides(result.rows), unmatchedRows };
  } catch {
    throw invalid();
  }
};

// Load replaces the edits of the committed request's records, as Load Label
// TSV replaces the label edits; edits of other records and of the other mode
// stay (R2).
export const replaceFeatureEdits = (featureOverrides, rows, mode, records) => {
  Object.values(featureOverrides).forEach((row) => {
    if (rowBelongsToRequest(row, mode, records)) updateFeatureOverride(featureOverrides, row, CLEARED_EDITS);
  });
  rows.forEach((row) => updateFeatureOverride(featureOverrides, { scope: mode, ...row }, editsOf(row)));
};

/**
 * @typedef {object} FeatureEditTableOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(value?: any) => { value: any }} ref Vue `ref`
 * @property {<T>(getter: () => T) => { value: T }} computed Vue `computed`
 * @property {() => Record<string, any> | null} getCommittedSession The committed canonical Session (render request and resources).
 * @property {((resourceId: string, kind: string) => any) | null} [readResourceRecordCount]
 *   Counts the records of a committed resource, for the Source recipe checks of the export.
 * @property {((payload: Record<string, any>) => Promise<{ result?: any }>) | null} [readFeatureOverrideTable]
 *   The diagram helper that reads a Feature Edits TSV (R7).
 * @property {() => any} projectFeatureEdits The root's projection of loaded feature edits onto the displayed Result (R3).
 */

/** @param {FeatureEditTableOptions} options */
export const createFeatureEditTableActions = ({
  state,
  ref,
  computed,
  getCommittedSession,
  readResourceRecordCount,
  readFeatureOverrideTable,
  projectFeatureEdits
}) => {
  const { featureOverrides } = state;

  const downloadFeatureEditTable = async () => {
    const committed = getCommittedSession();
    const records = committed?.renderRequest?.records || [];
    const rows = committed && records.length ? requestFeatureOverrides(featureOverrides, committed.renderRequest.mode, records) : [];
    if (!committed || rows.length === 0) {
      window.alert('No feature edits to export.');
      return false;
    }
    const { text, reason } = await buildFeatureOverrideTable({
      renderRequest: committed.renderRequest, resources: committed.resources, rows, readResourceRecordCount
    });
    if (!text) {
      window.alert(`Cannot export feature edits: ${reason}`);
      return false;
    }
    downloadTextFile('gbdraw_feature_override_table.tsv', text, 'text/tab-separated-values');
    const elsewhere = Object.values(featureOverrides).filter(hasEdit).length - rows.length;
    if (elsewhere > 0) {
      window.alert(`Exported ${rows.length} feature edit(s). ${elsewhere} feature edit(s) of records `
        + 'outside the current diagram were not exported.');
    }
    return true;
  };

  const failure = ref(null);
  const canRetryFeatureEditTableImport = computed(() => Boolean(failure.value
    && state.errorLog?.value === failure.value.error));
  const retryFeatureEditTableImport = () => (canRetryFeatureEditTableImport.value ? failure.value.retry() : false);
  const reselectFeatureEditTable = () => {
    if (canRetryFeatureEditTableImport.value) failure.value.input?.click?.();
  };

  const loadFeatureEditTable = async (event) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const input = event?.target;
    const file = input?.files?.[0];
    if (!file) return false;
    const committed = getCommittedSession();
    const records = committed?.renderRequest?.records || [];
    try {
      if (!committed || !records.length) {
        window.alert('No diagram is currently displayed. Generate a diagram, then load the feature edits TSV.');
        return false;
      }
      // app-setup wires the helper; a missing one throws here and is reported as a Feature edits TSV error.
      const response = await /** @type {NonNullable<typeof readFeatureOverrideTable>} */ (readFeatureOverrideTable)({
        files: [{ role: 'featureOverrides', bytes: await cloneFileBytesForTransfer(file) }],
        canonicalRequest: committed.renderRequest,
        resources: committed.resources
      });
      // A Generate that finished meanwhile names other records.
      if (getCommittedSession() !== committed || state.sessionOperationAvailability?.()) return false;
      const { rows, unmatchedRows } = admitFeatureOverrideTable(response?.result, records);
      replaceFeatureEdits(featureOverrides, rows, committed.renderRequest.mode, records);
      await projectFeatureEdits();
      if (state.errorLog?.value === failure.value?.error) state.errorLog.value = null;
      failure.value = null;
      let message = `Loaded ${rows.length + unmatchedRows.length} row(s). Applied ${rows.length} feature edit(s).`;
      if (unmatchedRows.length > 0) {
        message += ` ${unmatchedRows.length} row(s) name no record or feature of the current diagram and were not applied.`;
      }
      window.alert(message);
      return true;
    } catch (error) {
      const model = normalizeUserFacingError(error, { operation: OPERATION, stage: 'helper' });
      if (state.errorLog) state.errorLog.value = model;
      const sourceInput = event.sourceInput || input;
      failure.value = { error: model, input: sourceInput,
        retry: () => loadFeatureEditTable({ target: { files: [file], value: '' }, sourceInput }) };
      return false;
    } finally {
      if (input?.files?.[0] === file) input.value = '';
    }
  };

  return {
    canRetryFeatureEditTableImport,
    downloadFeatureEditTable,
    loadFeatureEditTable,
    reselectFeatureEditTable,
    retryFeatureEditTableImport
  };
};
