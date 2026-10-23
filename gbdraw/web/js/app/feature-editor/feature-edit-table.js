// @ts-check
/** @import { ShowNotice } from '../app-messages.js' */
import { buildFeatureOverrideTable } from '../run-info.js';
import {
  FEATURE_OVERRIDE_EDIT_FIELDS,
  canonicalFeatureOverrides,
  recordKeyBelongsToRequest,
  requestFeatureOverrides,
  rowBelongsToRequest,
  updateFeatureOverride
} from '../../services/feature-placement.js';
import { diagnosticError, normalizeUserFacingError } from '../../utils/error-normalization.js';
import { cloneFileBytesForTransfer } from '../../services/file-content-cache.js';
import { downloadTextFile } from '../../services/text-download.js';
import { featureEditsOf } from '../../services/drawing/feature-edits.js';
import { joinSentences, t, tn } from '../../utils/i18n.js';
/** @import { FeatureOverrideDraftRow } from '../../contracts/session.generated.js' */
/** @import { FeatureRequestRecord } from '../../services/feature-placement.js' */
/** @import { FeatureEditsSlice } from '../../services/drawing/feature-edits.js' */
/** @import { CanonicalFeatureOverride } from '../../contracts/request.generated.js' */
/** @import { LabelTableInput } from './label-actions.js' */

// Export and Load Feature Edits TSV (design Q4 6.4, Owner Q2 = A): the
// per-feature edits of the committed records as a --feature_override_table, so
// the CLI and the Python API read the same file. Python reads a loaded table
// (R4, R9); rows that name no record or feature of the committed records are
// counted and not applied (Owner Q3 = A).
const OPERATION = 'readFeatureOverrideTable';
const CLEARED_EDITS = Object.freeze(Object.fromEntries(FEATURE_OVERRIDE_EDIT_FIELDS.map((field) => [field, null])));
/** @param {Readonly<Record<string, unknown>>} row */
const editsOf = (row) => Object.fromEntries(FEATURE_OVERRIDE_EDIT_FIELDS.map((field) => [field, row[field]]));
/** @param {Readonly<Record<string, unknown>> | null | undefined} row */
const hasEdit = (row) => FEATURE_OVERRIDE_EDIT_FIELDS.some((field) => row?.[field] != null);

// The helper's rows must be request rows of the committed records.
/**
 * @param {{ rows?: unknown, unmatchedRows?: unknown } | null | undefined} result The helper's reply, checked below.
 * @param {FeatureRequestRecord[]} records
 */
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
// TSV replaces the label edits; edits of other records stay (R2).
/**
 * @param {FeatureEditsSlice} slice `drawing.slices.featureEdits`.
 * @param {Array<CanonicalFeatureOverride | FeatureOverrideDraftRow>} rows
 * @param {FeatureRequestRecord[]} records
 */
export const replaceFeatureEdits = (slice, rows, records) => {
  Object.values(featureEditsOf(slice).featureOverrides).forEach((row) => {
    if (rowBelongsToRequest(row, records)) updateFeatureOverride(slice, row, CLEARED_EDITS);
  });
  rows.forEach((row) => updateFeatureOverride(slice, row, editsOf(row)));
};

/**
 * @typedef {object} FeatureEditTableOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(value?: any) => { value: any }} ref Vue `ref`
 * @property {<T>(getter: () => T) => { value: T }} computed Vue `computed`
 * @property {() => Record<string, any> | null} getCommittedSession The committed canonical Session (render request and resources).
 * @property {((resourceId: string, kind: 'genbank' | 'fasta') => number | Promise<number>) | null} [readResourceRecordCount]
 *   Counts the records of a committed resource, for the Source recipe checks of the export.
 * @property {((payload: Record<string, any>) => Promise<{ result?: any }>) | null} [readFeatureOverrideTable]
 *   The diagram helper that reads a Feature Edits TSV (R7).
 * @property {() => any} projectFeatureEdits The root's projection of loaded feature edits onto the displayed Result (R3).
 * @property {ShowNotice} [showNotice] The root's notice banner: a table import or export reports there.
 */

/** @param {FeatureEditTableOptions} options */
export const createFeatureEditTableActions = ({
  state,
  ref,
  computed,
  getCommittedSession,
  readResourceRecordCount,
  readFeatureOverrideTable,
  projectFeatureEdits,
  showNotice = () => {}
}) => {

  const downloadFeatureEditTable = async () => {
    const drawing = state.activeDrawing();
    const committed = getCommittedSession();
    const records = committed?.renderRequest?.records || [];
    const rows = committed && records.length ? requestFeatureOverrides(drawing.featureOverrides, records) : [];
    if (!committed || rows.length === 0) {
      showNotice({ text: t('There are no feature edits to export.') });
      return false;
    }
    const { text, reason } = await buildFeatureOverrideTable({
      renderRequest: committed.renderRequest, resources: committed.resources, rows, readResourceRecordCount
    });
    if (!text) {
      showNotice({ text: t('Cannot export feature edits: {reason}', { reason: String(reason) }), tone: 'warning' });
      return false;
    }
    downloadTextFile('gbdraw_feature_override_table.tsv', text, 'text/tab-separated-values');
    const elsewhere = Object.values(drawing.featureOverrides).filter(hasEdit).length - rows.length;
    if (elsewhere > 0) {
      showNotice({ text: tn(rows.length, 'Exported {count} feature edit.', 'Exported {count} feature edits.'), items: [
        tn(elsewhere, '{count} feature edit of a record outside the current diagram was not exported.', '{count} feature edits of records outside the current diagram were not exported.')
      ] });
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

  /** @param {{ target: LabelTableInput | null, sourceInput?: LabelTableInput }} event A file input's change event, or the retry's stand-in for it. */
  const loadFeatureEditTable = async (event) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const input = event?.target;
    const file = input?.files?.[0];
    if (!file) return false;
    const committed = getCommittedSession();
    const records = committed?.renderRequest?.records || [];
    try {
      if (!committed || !records.length) {
        showNotice({ text: t('No diagram is shown. Generate a diagram, then load the Feature Edits TSV.'), tone: 'warning' });
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
      replaceFeatureEdits(drawing.slices.featureEdits, rows, records);
      await projectFeatureEdits();
      if (state.errorLog?.value === failure.value?.error) state.errorLog.value = null;
      failure.value = null;
      showNotice({
        text: joinSentences(tn(rows.length + unmatchedRows.length, 'Loaded {count} row.', 'Loaded {count} rows.'), tn(rows.length, 'Applied {count} feature edit.', 'Applied {count} feature edits.')),
        items: unmatchedRows.length > 0
          ? [tn(unmatchedRows.length, '{count} row was not applied: no matching record, or no matching feature in that record.', '{count} rows were not applied: no matching record, or no matching feature in that record.')]
          : []
      });
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
