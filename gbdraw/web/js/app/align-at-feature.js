// @ts-check
// Align at feature (Linear): Python resolves where `anchor_feature` places each
// record (gbdraw/web_support/feature_anchor.py): it rotates a complete circular
// record so that the feature starts at its left end, and moves the other records
// sideways. The page shows what happens to each record, and Apply writes Python's
// record translations and the rotated records' display starts once, as one
// History step, like Similarity alignment Apply. The page computes no start, no
// translation, and no reference record.
import { diagnosticError, normalizeUserFacingError } from '../utils/error-normalization.js';
import { msg, t } from '../utils/i18n.js';
import { canonicalRecordReverseComplement } from '../services/record-display-model.js';
import { recordKeyExpansion } from '../services/record-draw-selection.js';
import { computed, ref } from '../utils/reactivity.js';

/** @import { CommittedCanonicalSession, RecordPlacementStart } from '../services/session-request.js' */
/** @import { CanonicalRenderRequest } from '../contracts/request.generated.js' */
/** @import { HelperPayload } from '../services/diagram-generation.js' */
/** @import { Ref } from '../utils/reactivity.js' */

/**
 * What Align at feature does to one drawn record (`FeatureAnchorRecord` in
 * gbdraw/web_support/feature_anchor.py): `rotated` turns a complete circular
 * record so that the feature starts at its left end; `reference`, when no
 * record is rotated, keeps its place as the first record with the feature;
 * `aligned` moves so that the feature starts below the left end of the first
 * rotated record, or below the reference's feature; `unmatched` has no such
 * feature and keeps its place and start.
 * @typedef {'reference' | 'rotated' | 'aligned' | 'unmatched'} FeatureAnchorStatus
 * @typedef {object} FeatureAnchorRecord
 * @property {string} recordKey A drawn record's key: the request record's, or `<key>:<k>` for a member of a record
 *   that draws every record of its file.
 * @property {string} recordId
 * @property {FeatureAnchorStatus} status
 * @property {number | null} startCoordinate The record's display start in the resolved request: the new start of a
 *   rotated record, else its own, or null.
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
 * Draws a candidate of the committed Session and commits it as one History step under `label`
 * (`runCommittedCanonicalCandidate` in app/run-analysis.js, wired by the composition root).
 * `starts` are the display starts the candidate writes (Python's, for the records it rotated): the run also commits
 * them to the records' drafts, so the next Generate draws them, and one Undo restores the drafts with the Result.
 * @typedef {(run: { canonical: CommittedCanonicalSession, label: string, starts: readonly RecordPlacementStart[] }) => Promise<{ status?: string, error?: unknown } | null | undefined>} RecordPlacementRun
 */
/**
 * The candidate Session: `committed` with these record translations, orientations and display starts, and no
 * alignment plan (`projectCommittedSimilarityAlignment` in services/session-request.js).
 * @typedef {(input: { committed: CommittedCanonicalSession, plan: null, translations: FeatureAnchorTranslation[],
 *   orientations: { recordKey: string, reverseComplement: boolean }[], starts: readonly RecordPlacementStart[] }) => CommittedCanonicalSession} RecordPlacementProjection
 */
/**
 * @typedef {object} FeatureAnchorPorts
 * @property {() => CommittedCanonicalSession} materializedSession The committed Session with
 *   the active alignment plan materialized into its record translations, and no plan.
 * @property {() => boolean} hasAlignmentPlan
 * @property {(operation: string, payload: FeatureAnchorHelperPayload) => Promise<{ result?: unknown }>} runHelperOperation
 * @property {string} resolveOperation
 * @property {RecordPlacementRun} runPlacement Draws `canonical` and commits it as one History step under `label`.
 * @property {(reason: string) => void} notePlanReplaced
 * @property {RecordPlacementProjection} projectPlacement
 * @property {() => CanonicalRenderRequest | null} committedRequest The committed request (reactive): the notice after
 *   Apply holds only while its record translations and starts are the ones Apply wrote, so Undo retires it.
 * @property {() => number} sessionReplacements The Sessions the page has installed (reactive): a Load replaces the
 *   document, so it retires what Apply did, also when the loaded positions are the ones Apply wrote.
 */

const STATUSES = new Set(['reference', 'rotated', 'aligned', 'unmatched']);

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
/** @param {unknown} value @returns {value is number} */
const isFinite = (value) => typeof value === 'number' && Number.isFinite(value);
/** @param {unknown} value @returns {value is number} */
const isBase = (value) => typeof value === 'number' && Number.isSafeInteger(value) && value >= 1;

/**
 * The request record a drawn record's key belongs to: the record of that key, or, for a record that draws every
 * record of its file, one of its members `<key>:<k>`, which Python draws one by one.
 * @param {{ recordKey: string, cardinality: string }} record
 * @param {string} key
 */
const drawsRecordKey = (record, key) => key === record.recordKey
  || (record.cardinality === 'all' && recordKeyExpansion(key, record.recordKey) !== null);

/**
 * Check the helper reply against the request it was asked for: one entry per
 * drawn record, by record key, in the request's order, every request record
 * drawn (a record that draws every record of its file by its members, each
 * once), and the translations in the same order. A rotated record names its
 * new start.
 * @param {unknown} value
 * @param {readonly { recordKey: string, cardinality: string }[]} requestRecords
 * @returns {FeatureAnchorReply}
 */
export const validateFeatureAnchorReply = (value, requestRecords) => {
  const fail = () => { throw diagnosticError('RESULT_INVALID', {}, { stage: 'helper', operation: 'resolveFeatureAnchor' }); };
  if (!isObject(value) || typeof value.selector !== 'string'
    || !Array.isArray(value.records) || !Array.isArray(value.recordTranslations)) return fail();
  const records = value.records.map((/** @type {unknown} */ record) => {
    if (!isObject(record) || typeof record.recordKey !== 'string' || typeof record.recordId !== 'string'
      || !STATUSES.has(String(record.status)) || !isFinite(record.beforeX) || !isFinite(record.afterX)
      || !(record.status === 'rotated' ? isBase(record.startCoordinate)
        : record.startCoordinate === null || isBase(record.startCoordinate))) return fail();
    return /** @type {FeatureAnchorRecord} */ ({ recordKey: record.recordKey, recordId: record.recordId,
      status: record.status, startCoordinate: record.startCoordinate, beforeX: record.beforeX, afterX: record.afterX });
  });
  const keys = records.map(({ recordKey }) => recordKey);
  const owners = keys.map((key) => requestRecords.findIndex((record) => drawsRecordKey(record, key)));
  const members = keys.map((key, index) => (owners[index] < 0 ? null : recordKeyExpansion(key, requestRecords[owners[index]].recordKey)));
  // Each request record in its order, and the members of one in theirs.
  const ordered = owners.every((owner, index) => owner >= 0 && (index === 0 || owner > owners[index - 1]
    || (owner === owners[index - 1] && members[index] !== null && members[index - 1] !== null
      && Number(members[index]) > Number(members[index - 1]))));
  if (!ordered || new Set(owners).size !== requestRecords.length) return fail();
  const recordTranslations = value.recordTranslations.map((/** @type {unknown} */ entry, /** @type {number} */ index) => {
    if (!isObject(entry) || entry.recordKey !== keys[index] || !isFinite(entry.x) || !isFinite(entry.y)) return fail();
    return { recordKey: keys[index], x: entry.x, y: entry.y };
  });
  if (recordTranslations.length !== keys.length) return fail();
  return { selector: value.selector, records, recordTranslations };
};

/**
 * @param {FeatureAnchorRecord} record
 * @param {number} index
 */
const previewLine = (record, index) => {
  const name = `#${index + 1} ${record.recordId}`;
  if (record.status === 'rotated') return t('{record} starts at base {start} (rotated)', { record: name, start: String(record.startCoordinate) });
  if (record.status === 'reference') return t('{record} keeps its place (first record with the feature)', { record: name });
  if (record.status === 'unmatched') return t('{record} keeps its place (no matching feature)', { record: name });
  const distance = Math.round((record.afterX - record.beforeX) * 10) / 10;
  return distance === 0
    ? t('{record} is already aligned', { record: name })
    : t('{record} moves {distance} px', { record: name, distance: distance > 0 ? `+${distance}` : String(distance) });
};

/** @param {FeatureAnchorPorts} ports */
export const createFeatureAnchorAction = ({
  materializedSession, hasAlignmentPlan, runHelperOperation, resolveOperation, runPlacement, notePlanReplaced, projectPlacement,
  committedRequest, sessionReplacements
}) => {
  /** @type {Ref<string>} */
  const selector = ref('');
  /** @type {Ref<'idle' | 'resolving' | 'reviewing' | 'applying'>} */
  const status = ref('idle');
  /** @type {Ref<FeatureAnchorPreview | null>} */
  const preview = ref(null);
  /** @type {Ref<Record<string, unknown> | null>} */
  const error = ref(null);
  /** @type {Ref<{ text: string, written: FeatureAnchorTranslation[], starts: readonly RecordPlacementStart[], session: number }>} */
  const noticeState = ref({ text: '', written: [], starts: [], session: -1 });
  // Align at feature holds while the committed request draws the records where Apply put them, with the starts it
  // wrote, in the Session Apply wrote them to. Undo (also of an Apply that only rotated, whose positions may equal the
  // earlier ones), Reset record positions, a later placement, or a Session Load retires it; a drag with no plan
  // edits the composition, not the translations, so it does not. Its notice shows while it holds. Wrap Rows does not
  // read it: Align at feature counts there only through the records its Apply moved (PD-OI-113 scenario revisions
  // 2 and 3).
  const applied = computed(() => {
    const { written, starts, session } = noticeState.value;
    const request = committedRequest();
    const drawn = request && 'recordTranslations' in request.layout ? request.layout.recordTranslations : [];
    return request !== null && session === sessionReplacements() && written.length > 0 && drawn.length === written.length
      && written.every((entry, index) => {
        const current = drawn[index];
        return current.recordKey === entry.recordKey && current.x === entry.x && current.y === entry.y;
      })
      && starts.every(({ recordKey, startCoordinate }) => request.records.some((record) => record.recordKey === recordKey
        && record.display?.startCoordinate === startCoordinate));
  });
  const notice = computed(() => (applied.value ? noticeState.value.text : ''));
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
    noticeState.value = { text: '', written: [], starts: [], session: -1 };
    status.value = 'resolving';
    try {
      const canonical = materializedSession();
      const reply = await runHelperOperation(resolveOperation, {
        canonicalRequest: canonical.renderRequest, selector: selector.value,
        // The Worker stages the committed resource descriptors as Generate does.
        resources: /** @type {HelperPayload['resources']} */ (canonical.resources) });
      if (expected !== actionId) return { status: 'stale' };
      const checked = validateFeatureAnchorReply(reply?.result, canonical.renderRequest.records);
      // A record of the request that draws every record of its file is drawn as its members, one by one: Apply would
      // write each member's start and place, which one File card for all of them cannot hold (OV-2474).
      if (checked.records.some(({ recordKey }) => !canonical.renderRequest.records.some((record) => record.recordKey === recordKey))) {
        throw diagnosticError('INPUT_INVALID', { field: 'anchor_feature', reason: 'ANCHOR_RECORD_CARDS' });
      }
      preview.value = { reply: checked,
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

  // Write the previewed translations and the rotated records' starts as one History step, on the Result they were
  // resolved on.
  const apply = async () => {
    const current = preview.value;
    if (!current || status.value !== 'reviewing') return { status: 'rejected' };
    const expected = ++actionId;
    /** @type {CommittedCanonicalSession} */
    let canonical;
    /** @type {boolean} */
    let clearsAlignment;
    const starts = current.reply.records.flatMap(({ recordKey, status: recordStatus, startCoordinate }) => (
      recordStatus === 'rotated' && startCoordinate !== null ? [{ recordKey, startCoordinate }] : []));
    // A Result that changed or is gone since the preview (a live edit, a Load) refuses Apply
    // with the reason, and the action returns to idle.
    try {
      /** @type {CommittedCanonicalSession | null} */
      let committed = null;
      try { committed = materializedSession(); } catch { committed = null; }
      if (!committed || JSON.stringify(committed.renderRequest) !== current.request) {
        throw diagnosticError('INPUT_INVALID', { field: 'anchor_feature', reason: 'ANCHOR_STALE' });
      }
      clearsAlignment = hasAlignmentPlan();
      canonical = projectPlacement({ committed, plan: null,
        translations: current.reply.recordTranslations,
        orientations: committed.renderRequest.records.map((record) => ({ recordKey: record.recordKey,
          reverseComplement: canonicalRecordReverseComplement(record) })),
        starts });
    } catch (cause) {
      cancel();
      publishError(cause);
      return { status: 'refused' };
    }
    status.value = 'applying';
    /** @type {{ status?: string, error?: unknown } | null | undefined} */
    let outcome;
    try { outcome = await runPlacement({ canonical, label: msg('Align at feature'), starts }); }
    catch (cause) { outcome = { status: 'error', error: cause }; }
    if (expected !== actionId) return { status: 'stale' };
    if (outcome?.status === 'ok') {
      if (clearsAlignment) notePlanReplaced(msg('records aligned at a feature.'));
      noticeState.value = { text: t('Records aligned at {selector}.', { selector: current.reply.selector }),
        written: current.reply.recordTranslations, starts, session: sessionReplacements() };
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
