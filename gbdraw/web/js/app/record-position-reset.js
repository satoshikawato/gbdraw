// @ts-check
// Reset record positions (Linear, PD-OI-113 scenario revision 2): one button under a Wrap Rows field that a drawn
// record moved sideways blocks, offered only while no comparison blocks it (`offersRecordPositionReset`, D-74). It
// writes every record's sideways position (x) as 0 through the channel Align at feature's Apply writes them through:
// a candidate of the committed Session with those record translations, drawn and committed as one History step
// (`runPlacement`), so one Undo restores the positions. A vertical position (y) and the record directions stay.
// Align at feature's hold follows the positions: the press retires it and its Undo brings it back.
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { msg } from '../utils/i18n.js';
import { canonicalRecordReverseComplement } from '../services/record-display-model.js';
import { recordTranslationsFor } from '../services/linear-record-layout.js';
import { ref } from '../utils/reactivity.js';

/** @import { CommittedCanonicalSession } from '../services/session-request.js' */
/** @import { CanonicalRenderRequest } from '../contracts/request.generated.js' */
/** @import { FeatureAnchorTranslation, RecordPlacementRun } from './align-at-feature.js' */

/**
 * @typedef {object} RecordPositionResetPorts
 * @property {() => CommittedCanonicalSession} materializedSession The committed Session (Align at feature's port).
 * @property {RecordPlacementRun} runPlacement Draws `canonical` and commits it as one History step under `label`
 *   (Align at feature's port).
 * @property {(input: { committed: CommittedCanonicalSession, plan: null, translations: FeatureAnchorTranslation[],
 *   orientations: { recordKey: string, reverseComplement: boolean }[] }) => CommittedCanonicalSession} projectPlacement
 * @property {(error: Record<string, unknown>) => void} onError Shows why the positions were not reset.
 */

/**
 * Every record of the committed request at x 0, its y kept: one entry per record, by record key, read from the
 * translations the request carries (`recordTranslationsFor`).
 * @param {CanonicalRenderRequest} request
 * @returns {FeatureAnchorTranslation[]}
 */
export const recordTranslationsAtX0 = (request) => {
  const stored = 'recordTranslations' in request.layout ? request.layout.recordTranslations : [];
  return request.records.map(({ recordKey }) => {
    const [entry] = recordTranslationsFor(stored, [recordKey]);
    return { recordKey, x: 0, y: typeof entry?.y === 'number' ? entry.y : 0 };
  });
};

/** @param {RecordPositionResetPorts} ports */
export const createRecordPositionReset = ({ materializedSession, runPlacement, projectPlacement, onError }) => {
  const busy = ref(false);
  /** @param {unknown} cause */
  const publishError = (cause) => onError(/** @type {Record<string, unknown>} */ (
    normalizeUserFacingError(cause || { code: 'UNKNOWN' }, { operation: 'generate', stage: 'helper' })));
  const reset = async () => {
    if (busy.value) return { status: 'busy' };
    /** @type {CommittedCanonicalSession} */
    let canonical;
    try {
      const committed = materializedSession();
      canonical = projectPlacement({ committed, plan: null,
        translations: recordTranslationsAtX0(committed.renderRequest),
        orientations: committed.renderRequest.records.map((record) => ({ recordKey: record.recordKey,
          reverseComplement: canonicalRecordReverseComplement(record) })) });
    } catch (cause) {
      publishError(cause);
      return { status: 'refused' };
    }
    busy.value = true;
    try {
      /** @type {{ status?: string, error?: unknown } | null | undefined} */
      let outcome;
      try { outcome = await runPlacement({ canonical, label: msg('Reset record positions') }); }
      catch (cause) { outcome = { status: 'error', error: cause }; }
      if (outcome?.status === 'error') publishError(outcome.error);
      return { status: String(outcome?.status || 'error') };
    } finally {
      busy.value = false;
    }
  };
  return { busy, reset };
};
