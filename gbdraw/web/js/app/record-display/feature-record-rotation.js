import { resolveFeatureAnchor } from './feature-anchor.js';

/**
 * Coordinate one explicit popup feature through the existing domain owners.
 * The controller owns no request editing, generation, admission, or History stack.
 */
export const createFeatureRecordRotationAction = ({
  recordDisplayControls,
  getCommittedSession,
  projectCommittedRecordTransform,
  runCommittedCanonicalCandidate
}) => {
  if (!recordDisplayControls
    || typeof getCommittedSession !== 'function'
    || typeof projectCommittedRecordTransform !== 'function'
    || typeof runCommittedCanonicalCandidate !== 'function') {
    throw new Error('Feature record rotation owners are unavailable.');
  }

  const resolve = ({ feature, intent }) => {
    const { row, target } = recordDisplayControls.targetForFeature(feature);
    const resolved = resolveFeatureAnchor({
      recordLength: target.recordLength,
      effectiveCircular: target.effectiveCircular,
      cropped: target.cropped,
      currentReverseComplement: target.committedReverseComplement,
      identity: {
        recordKey: target.recordKey,
        biologicalFeatureId: feature?.biological_feature_id
      },
      parts: feature?.location_parts,
      profile: feature?.anchorProfile,
      intent
    });
    return { row, target, resolved };
  };

  const apply = async ({ feature, intent }) => {
    const { row, target, resolved } = resolve({ feature, intent });
    if (!resolved.eligibility.enabled) {
      throw new Error(resolved.eligibility.message);
    }
    const projection = projectCommittedRecordTransform({
      committed: getCommittedSession(),
      target,
      transform: {
        recordLength: target.recordLength,
        startCoordinate: resolved.startCoordinate,
        reverseComplement: resolved.reverseComplement
      }
    });
    const outcome = await runCommittedCanonicalCandidate({
      canonical: projection.canonical,
      captureIntentCheckpoint: () => recordDisplayControls.captureTargetDraft(row),
      restoreIntentCheckpoint: (checkpoint) => (
        recordDisplayControls.restoreTargetDraft(checkpoint)
      ),
      commitIntent: () => recordDisplayControls.commitResolvedTransform(row, {
        startCoordinate: resolved.startCoordinate,
        reverseComplement: resolved.reverseComplement,
        anchorIntent: resolved.provenance
      })
    });
    return { ...outcome, receipt: projection.receipt, resolved };
  };

  return Object.freeze({ resolve, apply });
};
