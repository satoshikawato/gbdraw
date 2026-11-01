// @ts-check
// Align at feature (Linear): Python resolves where `anchor_feature` places each
// record (gbdraw/web_support/feature_anchor.py), the page shows which records
// move, and Apply writes the record translations once, as one History step,
// like Similarity alignment Apply. The record translations stay the only owner
// of a record's place; the page copies no start or translation rule.
import { diagnosticError, normalizeUserFacingError } from '../utils/error-normalization.js';
import { msg, t } from '../utils/i18n.js';
import { canonicalRecordReverseComplement } from '../services/record-display-model.js';
import { computed, ref } from '../utils/reactivity.js';

/** @import { CommittedCanonicalSession } from '../services/session-request.js' */
/** @import { HelperPayload } from '../services/diagram-generation.js' */
/** @import { Ref } from '../utils/reactivity.js' */

/**
 * What Align at feature does to one drawn record (`FeatureAnchorRecord` in
 * gbdraw/web_support/feature_anchor.py): `reference` keeps its place as the
 * first record with the feature, `aligned` moves so that the feature starts
 * below the reference's, `unmatched` has no such feature and keeps its place.
 * @typedef {'reference' | 'aligned' | 'unmatched'} FeatureAnchorStatus
 * @typedef {object} FeatureAnchorRecord
 * @property {string} recordKey
 * @property {string} recordId
 * @property {FeatureAnchorStatus} status
 * @property {number} beforeX
 * @property {number} afterX
 */
/**
 * One `renderRequest.layout.recordTranslations` entry (`EncodedRecordTranslation`).
 * @typedef {{ recordKey: string, x: number, y: number }} FeatureAnchorTranslation
 */
/**
 * The `resolveFeatureAnchor` reply (`FeatureAnchorReply`).
 * @typedef {object} FeatureAnchorReply
 * @property {string} selector
 * @property {FeatureAnchorRecord[]} records
 * @property {FeatureAnchorTranslation[]} recordTranslations
 */
/**
 * The `resolveFeatureAnchor` request: the committed request, the selector as typed, and the
 * committed resources the page stages for the Worker.
 * @typedef {object} FeatureAnchorHelperPayload
 * @property {CommittedCanonicalSession['renderRequest']} canonicalRequest
 * @property {string} selector
 * @property {HelperPayload['resources']} resources
 */
/**
 * @typedef {object} FeatureAnchorPreview A resolved selector, waiting for Apply or Cancel.
 * @property {FeatureAnchorReply} reply
 * @property {string} request The request it was resolved on, as JSON: Apply refuses another.
 * @property {boolean} clearsAlignment Whether Apply clears the active Similarity alignment.
 */
/**
 * @typedef {object} FeatureAnchorPorts
 * @property {() => CommittedCanonicalSession} materializedSession The committed Session with
 *   the active alignment plan materialized into its record translations, and no plan.
 * @property {() => boolean} hasAlignmentPlan
 * @property {(operation: string, payload: FeatureAnchorHelperPayload) => Promise<{ result?: unknown }>} runHelperOperation
 * @property {string} resolveOperation
 * @property {(run: { canonical: CommittedCanonicalSession, label: string }) => Promise<{ status?: string, error?: unknown } | null | undefined>} runPlacement
 *   Draws `canonical` and commits it as one History step under `label`.
 * @property {(reason: string) => void} notePlanReplaced
 * @property {(input: { committed: CommittedCanonicalSession, plan: null, translations: FeatureAnchorTranslation[],
 *   orientations: { recordKey: string, reverseComplement: boolean }[] }) => CommittedCanonicalSession} projectPlacement
 *   The candidate Session: `committed` with these record translations and orientations, and no alignment plan.
 */

const STATUSES = new Set(['reference', 'aligned', 'unmatched']);

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
/** @param {unknown} value @returns {value is number} */
const isFinite = (value) => typeof value === 'number' && Number.isFinite(value);

/**
 * Check the helper reply against the request it was asked for: one entry per
 * record, in the request's order, by record key.
 * @param {unknown} value
 * @param {readonly string[]} recordKeys
 * @returns {FeatureAnchorReply}
 */
export const validateFeatureAnchorReply = (value, recordKeys) => {
  const fail = () => { throw diagnosticError('RESULT_INVALID', {}, { stage: 'helper', operation: 'resolveFeatureAnchor' }); };
  if (!isObject(value) || typeof value.selector !== 'string'
    || !Array.isArray(value.records) || !Array.isArray(value.recordTranslations)) return fail();
  const records = value.records.map((/** @type {unknown} */ record, /** @type {number} */ index) => {
    if (!isObject(record) || record.recordKey !== recordKeys[index] || typeof record.recordId !== 'string'
      || !STATUSES.has(String(record.status)) || !isFinite(record.beforeX) || !isFinite(record.afterX)) return fail();
    return /** @type {FeatureAnchorRecord} */ ({ recordKey: record.recordKey, recordId: record.recordId,
      status: record.status, beforeX: record.beforeX, afterX: record.afterX });
  });
  const recordTranslations = value.recordTranslations.map((/** @type {unknown} */ entry, /** @type {number} */ index) => {
    if (!isObject(entry) || entry.recordKey !== recordKeys[index] || !isFinite(entry.x) || !isFinite(entry.y)) return fail();
    return { recordKey: recordKeys[index], x: entry.x, y: entry.y };
  });
  if (records.length !== recordKeys.length || recordTranslations.length !== recordKeys.length) return fail();
  return { selector: value.selector, records, recordTranslations };
};

/**
 * @param {FeatureAnchorRecord} record
 * @param {number} index
 */
const previewLine = (record, index) => {
  const name = `#${index + 1} ${record.recordId}`;
  if (record.status === 'reference') return t('{record}: keeps its place (first record with the feature)', { record: name });
  if (record.status === 'unmatched') return t('{record}: keeps its place (no matching feature)', { record: name });
  const distance = Math.round((record.afterX - record.beforeX) * 10) / 10;
  return distance === 0
    ? t('{record}: already aligned', { record: name })
    : t('{record}: moves {distance} px', { record: name, distance: distance > 0 ? `+${distance}` : String(distance) });
};

/** @param {FeatureAnchorPorts} ports */
export const createFeatureAnchorAction = ({
  materializedSession, hasAlignmentPlan, runHelperOperation, resolveOperation, runPlacement, notePlanReplaced, projectPlacement
}) => {
  /** @type {Ref<string>} */
  const selector = ref('');
  /** @type {Ref<'idle' | 'resolving' | 'reviewing' | 'applying'>} */
  const status = ref('idle');
  /** @type {Ref<FeatureAnchorPreview | null>} */
  const preview = ref(null);
  /** @type {Ref<Record<string, unknown> | null>} */
  const error = ref(null);
  /** @type {Ref<string>} */
  const notice = ref('');
  let actionId = 0;
  const busy = computed(() => status.value === 'resolving' || status.value === 'applying');
  const previewLines = computed(() => (preview.value ? preview.value.reply.records.map(previewLine) : []));

  /** @param {unknown} cause */
  const publishError = (cause) => {
    error.value = /** @type {Record<string, unknown>} */ (
      normalizeUserFacingError(cause || { code: 'UNKNOWN' }, { operation: resolveOperation, stage: 'helper' }));
  };

  const cancel = () => {
    actionId += 1;
    preview.value = null;
    status.value = 'idle';
  };

  // Resolve the selector on the Result as it is drawn; nothing changes until Apply.
  const resolve = async () => {
    if (busy.value) return { status: 'busy' };
    const expected = ++actionId;
    preview.value = null;
    error.value = null;
    notice.value = '';
    status.value = 'resolving';
    try {
      const canonical = materializedSession();
      const reply = await runHelperOperation(resolveOperation, {
        canonicalRequest: canonical.renderRequest, selector: selector.value,
        // The Worker stages the committed resource descriptors as Generate does.
        resources: /** @type {HelperPayload['resources']} */ (canonical.resources) });
      if (expected !== actionId) return { status: 'stale' };
      const recordKeys = canonical.renderRequest.records.map(({ recordKey }) => recordKey);
      preview.value = { reply: validateFeatureAnchorReply(reply?.result, recordKeys),
        request: JSON.stringify(canonical.renderRequest), clearsAlignment: hasAlignmentPlan() };
      status.value = 'reviewing';
      return { status: 'reviewing' };
    } catch (cause) {
      if (expected !== actionId) return { status: 'stale' };
      status.value = 'idle';
      publishError(cause);
      return { status: 'error' };
    }
  };

  // Write the previewed translations as one History step, on the Result they were resolved on.
  const apply = async () => {
    const current = preview.value;
    if (!current || status.value !== 'reviewing') return { status: 'rejected' };
    const expected = ++actionId;
    const committed = materializedSession();
    if (JSON.stringify(committed.renderRequest) !== current.request) {
      cancel();
      publishError(diagnosticError('INPUT_INVALID', { field: 'anchor_feature', reason: 'ANCHOR_STALE' }));
      return { status: 'stale' };
    }
    const clearsAlignment = hasAlignmentPlan();
    const canonical = projectPlacement({ committed, plan: null,
      translations: current.reply.recordTranslations,
      orientations: committed.renderRequest.records.map((record) => ({ recordKey: record.recordKey,
        reverseComplement: canonicalRecordReverseComplement(record) })) });
    status.value = 'applying';
    /** @type {{ status?: string, error?: unknown } | null | undefined} */
    let outcome;
    try { outcome = await runPlacement({ canonical, label: msg('Align at feature') }); }
    catch (cause) { outcome = { status: 'error', error: cause }; }
    if (expected !== actionId) return { status: 'stale' };
    if (outcome?.status === 'ok') {
      if (clearsAlignment) notePlanReplaced(msg('records aligned at a feature.'));
      notice.value = t('Records aligned at {selector}.', { selector: current.reply.selector });
      preview.value = null;
      status.value = 'idle';
      return { status: 'ok' };
    }
    status.value = 'reviewing';
    if (outcome?.status === 'error') publishError(outcome.error);
    return { status: String(outcome?.status || 'error') };
  };

  return { selector, status, busy, preview, previewLines, error, notice, resolve, apply, cancel };
};
