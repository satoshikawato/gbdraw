// @ts-check
// Reset record positions (Linear, PD-OI-113 scenario revision 2, Q5 A): one button under a Wrap Rows field that a
// drawn record moved sideways, or a Similarity alignment, blocks (`offersRecordPositionReset`). It writes every
// record's sideways position (x) as 0 through the channel Align at feature's Apply writes them through: a candidate
// of the committed Session (the alignment plan's positions materialized) with those record translations and no
// plan, drawn and committed as one History step (`runPlacement`), so one Undo restores the positions and the plan.
// A vertical position (y), the record directions, and the display starts (rotations, PD-OI-113 scenario revision 3)
// stay.
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { msg } from '../utils/i18n.js';
import { canonicalRecordReverseComplement } from '../services/record-display-model.js';
import { recordTranslationsFor } from '../services/linear-record-layout.js';
import { ref } from '../utils/reactivity.js';

/** @import { CommittedCanonicalSession } from '../services/session-request.js' */
/** @import { UserFacingError } from '../utils/error-normalization.js' */
/** @import { CanonicalRenderRequest } from '../contracts/request.generated.js' */
/** @import { FeatureAnchorTranslation, RecordPlacementProjection, RecordPlacementRun } from './align-at-feature.js' */

/**
 * @typedef {object} RecordPositionResetPorts
 * @property {() => CommittedCanonicalSession} materializedSession The committed Session with the active alignment
 *   plan materialized into its record translations, and no plan; from a Result an older gbdraw drew, the committed
 *   translations (the Similarity alignment's `repairBaseSession`).
 * @property {() => boolean} hasAlignmentPlan
 * @property {(reason: string) => void} notePlanReplaced Says that the press ended the Similarity alignment.
 * @property {RecordPlacementRun} runPlacement Draws `canonical` and commits it as one History step under `label`
 *   (Align at feature's port).
 * @property {RecordPlacementProjection} projectPlacement
 * @property {(error: UserFacingError) => void} onError Shows why the positions were not reset.
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
export const createRecordPositionReset = ({
  materializedSession, hasAlignmentPlan, notePlanReplaced, runPlacement, projectPlacement, onError
}) => {
  const busy = ref(false);
  /** @param {unknown} cause */
  // The value is always truthy, and only a falsy value normalizes to null.
  const publishError = (cause) => onError(/** @type {UserFacingError} */ (
    normalizeUserFacingError(cause || { code: 'UNKNOWN' }, { operation: 'generate', stage: 'helper' })));
  const reset = async () => {
    if (busy.value) return { status: 'busy' };
    /** @type {CommittedCanonicalSession} */
    let canonical;
    let endsAlignment = false;
    try {
      const committed = materializedSession();
      endsAlignment = hasAlignmentPlan();
      canonical = projectPlacement({ committed, plan: null,
        translations: recordTranslationsAtX0(committed.renderRequest),
        orientations: committed.renderRequest.records.map((record) => ({ recordKey: record.recordKey,
          reverseComplement: canonicalRecordReverseComplement(record) })),
        // Sideways only: the display starts stay, so a record Align at feature rotated keeps its rotation.
        starts: [] });
    } catch (cause) {
      publishError(cause);
      return { status: 'refused' };
    }
    busy.value = true;
    try {
      /** @type {{ status?: string, error?: unknown } | null | undefined} */
      let outcome;
      try { outcome = await runPlacement({ canonical, label: msg('Reset record positions'), starts: [] }); }
      catch (cause) { outcome = { status: 'error', error: cause }; }
      if (outcome?.status === 'ok' && endsAlignment) notePlanReplaced(msg('record positions reset.'));
      else if (outcome?.status === 'error') publishError(outcome.error);
      return { status: String(outcome?.status || 'error') };
    } finally {
      busy.value = false;
    }
  };
  return { busy, reset };
};
