import { resolveFeatureAnchor } from './feature-anchor.js';

const featureIdentity = (feature) => ({
  recordKey: String(feature?.record_key ?? feature?.recordKey ?? ''),
  biologicalFeatureId: String(
    feature?.biological_feature_id ?? feature?.biologicalFeatureId ?? ''
  )
});

const sameFeatureIdentity = (feature, identity) => {
  const candidate = featureIdentity(feature);
  return candidate.recordKey === identity.recordKey
    && candidate.biologicalFeatureId === identity.biologicalFeatureId;
};

/**
 * Coordinate one explicit popup feature through the existing domain owners.
 * The controller owns no request editing, generation, admission, or History stack.
 */
export const createFeatureRecordRotationAction = ({
  recordDisplayControls,
  getCommittedSession,
  projectCommittedRecordTransform,
  runCommittedCanonicalCandidate,
  resolveCurrentFeature = null,
  isCurrentFeature = null
}) => {
  if (!recordDisplayControls
    || typeof getCommittedSession !== 'function'
    || typeof projectCommittedRecordTransform !== 'function'
    || typeof runCommittedCanonicalCandidate !== 'function') {
    throw new Error('Feature record rotation owners are unavailable.');
  }

  const currentFeature = (feature) => {
    const identity = featureIdentity(feature);
    if (!identity.recordKey || !identity.biologicalFeatureId) {
      throw new Error('Feature rotation requires a stable record and feature identity.');
    }
    const resolved = typeof resolveCurrentFeature === 'function'
      ? resolveCurrentFeature(identity)
      : feature;
    if (!resolved || !sameFeatureIdentity(resolved, identity)) {
      throw new Error('The popup feature is no longer present in the current Result.');
    }
    if (typeof isCurrentFeature === 'function' && !isCurrentFeature(resolved)) {
      throw new Error('The popup feature source changed after the popup opened.');
    }
    return resolved;
  };

  const resolve = ({ feature, intent }) => {
    const resolvedFeature = currentFeature(feature);
    const { row, target } = recordDisplayControls.targetForFeature(resolvedFeature);
    const resolved = resolveFeatureAnchor({
      recordLength: target.recordLength,
      effectiveCircular: target.effectiveCircular,
      cropped: target.cropped,
      currentReverseComplement: target.committedReverseComplement,
      identity: {
        recordKey: target.recordKey,
        biologicalFeatureId: featureIdentity(resolvedFeature).biologicalFeatureId
      },
      parts: resolvedFeature?.location_parts,
      profile: resolvedFeature?.anchorProfile,
      intent
    });
    return { feature: resolvedFeature, row, target, resolved };
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

  return Object.freeze({ currentFeature, resolve, apply });
};

const parseSignedOffset = (value) => {
  const text = String(value ?? '').trim();
  if (!text) {
    return { valid: false, value: Number.NaN, message: 'Enter a signed whole-number offset.' };
  }
  if (!/^[+-]?\d+$/.test(text)) {
    return {
      valid: false,
      value: Number.NaN,
      message: 'Offset must be a signed whole number of base pairs.'
    };
  }
  const number = Number(text);
  if (!Number.isSafeInteger(number)) {
    return {
      valid: false,
      value: Number.NaN,
      message: 'Offset is outside the supported integer range.'
    };
  }
  return { valid: true, value: number, message: '' };
};

const emptyCapability = (message = '') => ({ enabled: false, code: null, message });

const initialDraft = () => ({
  active: false,
  feature: null,
  featureLabel: '',
  identity: { recordKey: '', biologicalFeatureId: '' },
  recordLabel: '',
  placement: 'anchor',
  anchor: 'five-prime',
  offsetText: '0',
  orientForward: false,
  exactFeatureEnd: false,
  anchorChoices: [
    { value: 'five-prime', label: '5′ end', ...emptyCapability() },
    { value: 'midpoint', label: 'Midpoint', ...emptyCapability() },
    { value: 'three-prime', label: '3′ end', ...emptyCapability() }
  ],
  orientCapability: emptyCapability(),
  featureEndCapability: emptyCapability(),
  offsetError: '',
  disabledReason: '',
  canApply: false,
  pending: false,
  status: '',
  statusKind: 'idle',
  startCoordinate: null,
  orientationLabel: 'unchanged',
  displayedStrandLabel: '',
  placementLabel: 'Custom anchor'
});

const copyCapability = (capability) => ({
  enabled: Boolean(capability?.enabled),
  code: capability?.code || null,
  message: String(capability?.message || '')
});

/**
 * Own only the popup's ephemeral draft and presentation state. Domain math and
 * transactional generation stay in createFeatureRecordRotationAction.
 */
export const createFeatureRecordRotationWorkflow = ({
  action,
  makeReactive = (value) => value,
  onRebind = null
}) => {
  if (!action || typeof action.resolve !== 'function' || typeof action.apply !== 'function') {
    throw new Error('Feature record rotation action is unavailable.');
  }
  const draft = makeReactive(initialDraft());

  const intent = () => {
    const offset = parseSignedOffset(draft.offsetText);
    return {
      placement: draft.placement,
      anchor: draft.placement === 'anchor' ? draft.anchor : null,
      offsetBp: offset.value,
      orientForward: Boolean(draft.orientForward)
    };
  };

  const recompute = () => {
    if (!draft.active || !draft.feature) return null;
    const offset = parseSignedOffset(draft.offsetText);
    draft.offsetError = offset.message;
    try {
      const snapshot = action.resolve({ feature: draft.feature, intent: intent() });
      draft.feature = snapshot.feature;
      draft.identity = featureIdentity(snapshot.feature);
      draft.recordLabel = String(snapshot.target.recordId || snapshot.target.recordKey);
      draft.anchorChoices = [
        ['five-prime', '5′ end'],
        ['midpoint', 'Midpoint'],
        ['three-prime', '3′ end']
      ].map(([value, label]) => ({
        value,
        label,
        ...copyCapability(snapshot.resolved.capabilities.anchors[value])
      }));
      draft.orientCapability = copyCapability(snapshot.resolved.capabilities.orientForward);
      draft.featureEndCapability = copyCapability(snapshot.resolved.capabilities.featureEnd);
      draft.disabledReason = offset.message || snapshot.resolved.eligibility.message || '';
      draft.canApply = offset.valid && snapshot.resolved.eligibility.enabled && !draft.pending;
      draft.startCoordinate = snapshot.resolved.startCoordinate;
      draft.orientationLabel = snapshot.resolved.reverseComplement
        === snapshot.target.committedReverseComplement
        ? 'unchanged'
        : snapshot.resolved.reverseComplement
          ? 'reverse-complemented'
          : 'forward';
      const before = snapshot.resolved.displayedStrand.before;
      const after = snapshot.resolved.displayedStrand.after;
      draft.displayedStrandLabel = ['+', '-'].includes(before) && ['+', '-'].includes(after)
        ? `${before} → ${after}`
        : '';
      draft.placementLabel = draft.placement === 'feature-end'
        ? draft.exactFeatureEnd
          ? 'Exactly at feature end'
          : 'Feature end with custom offset'
        : 'Custom anchor';
      return snapshot;
    } catch (error) {
      draft.anchorChoices = draft.anchorChoices.map((choice) => ({
        ...choice,
        ...emptyCapability(error.message)
      }));
      draft.orientCapability = emptyCapability(error.message);
      draft.featureEndCapability = emptyCapability(error.message);
      draft.disabledReason = offset.message || error.message;
      draft.canApply = false;
      draft.startCoordinate = null;
      draft.orientationLabel = 'unchanged';
      draft.displayedStrandLabel = '';
      return null;
    }
  };

  const open = ({ feature, featureLabel = '' }) => {
    Object.assign(draft, initialDraft(), {
      active: true,
      feature,
      featureLabel: String(featureLabel || feature?.label || feature?.type || ''),
      identity: featureIdentity(feature)
    });
    recompute();
    return draft;
  };

  const close = () => {
    Object.assign(draft, initialDraft());
  };

  const setAnchor = (anchor) => {
    draft.placement = 'anchor';
    draft.anchor = String(anchor || '');
    draft.exactFeatureEnd = false;
    recompute();
  };

  const setOffset = (value) => {
    draft.offsetText = String(value ?? '');
    if (draft.placement === 'feature-end') draft.exactFeatureEnd = false;
    recompute();
  };

  const setOrientForward = (value) => {
    draft.orientForward = Boolean(value);
    recompute();
  };

  const placeAtFeatureEnd = () => {
    draft.placement = 'feature-end';
    draft.anchor = null;
    draft.offsetText = '0';
    draft.exactFeatureEnd = true;
    recompute();
  };

  const apply = async () => {
    if (draft.pending) return { status: 'pending' };
    const snapshot = recompute();
    if (!snapshot || !draft.canApply) {
      draft.status = draft.disabledReason || 'Record rotation is unavailable.';
      draft.statusKind = 'error';
      return { status: 'disabled' };
    }
    draft.pending = true;
    draft.canApply = false;
    draft.status = 'Regenerating the target record…';
    draft.statusKind = 'pending';
    try {
      const outcome = await action.apply({
        feature: snapshot.feature,
        intent: intent()
      });
      if (outcome?.status === 'ok') {
        const rebound = action.currentFeature(snapshot.feature);
        draft.feature = rebound;
        if (typeof onRebind === 'function') onRebind(rebound, draft.identity);
        draft.status = 'Record rotation applied and regenerated.';
        draft.statusKind = 'success';
      } else if (outcome?.status === 'canceled') {
        draft.status = 'Record rotation canceled. The previous Result was kept.';
        draft.statusKind = 'idle';
      } else {
        draft.status = 'Record rotation failed. The previous Result was kept.';
        draft.statusKind = 'error';
      }
      return outcome;
    } catch (error) {
      draft.disabledReason = error.message;
      draft.status = `${error.message} The previous Result was kept.`;
      draft.statusKind = 'error';
      return { status: 'error', error };
    } finally {
      draft.pending = false;
      recompute();
    }
  };

  return Object.freeze({
    draft,
    open,
    close,
    cancel: close,
    recompute,
    setAnchor,
    setOffset,
    setOrientForward,
    placeAtFeatureEnd,
    apply
  });
};
