import { resolveFeatureAnchor } from './feature-anchor.js';
import { RECORD_TARGET_NOT_DISCOVERED } from '../record-display-options.js';

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
  isCurrentFeature = null,
  readRecords = null
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

  // The explicit record action that reads records a Session Load left unread.
  const readTargetRecords = async () => (typeof readRecords === 'function'
    ? readRecords()
    : { status: 'unavailable', reason: 'Records cannot be read here.' });

  return Object.freeze({ currentFeature, resolve, apply, readRecords: readTargetRecords });
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
  // '' when the popup feature resolves to one record target that a control can
  // act on; otherwise the kind of the whole-target failure, whose one reason is
  // disabledReason and which hides the form.
  targetFailure: '',
  disabledReason: '',
  canApply: false,
  pending: false,
  reading: false,
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
  // Bumped whenever the draft is reset, so a late record read cannot write
  // into a closed or retargeted popup.
  let draftGeneration = 0;

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
      const { anchors, orientForward, featureEnd } = snapshot.resolved.capabilities;
      const controls = [...Object.values(anchors), orientForward, featureEnd];
      // A limit that disables every control (for example, a linear record) is a
      // whole-target limit: the reason line states it once (F4).
      const targetLimit = controls.every((capability) => !capability?.enabled
        && capability?.code === controls[0]?.code) ? controls[0] : null;
      draft.targetFailure = targetLimit ? String(targetLimit.code || 'unavailable') : '';
      const controlCapability = (capability) => ({
        ...copyCapability(capability),
        ...(targetLimit ? { message: '' } : {})
      });
      draft.anchorChoices = [
        ['five-prime', '5′ end'],
        ['midpoint', 'Midpoint'],
        ['three-prime', '3′ end']
      ].map(([value, label]) => ({
        value,
        label,
        ...controlCapability(anchors[value])
      }));
      draft.orientCapability = controlCapability(orientForward);
      draft.featureEndCapability = controlCapability(featureEnd);
      draft.disabledReason = targetLimit?.message
        || offset.message || snapshot.resolved.eligibility.message || '';
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
      // A whole-target failure has one reason; per-control messages are only
      // for control-specific limits of a resolved target.
      draft.targetFailure = error?.kind || 'unavailable';
      draft.anchorChoices = draft.anchorChoices.map((choice) => ({
        ...choice,
        ...emptyCapability()
      }));
      draft.orientCapability = emptyCapability();
      draft.featureEndCapability = emptyCapability();
      draft.disabledReason = error.message;
      draft.canApply = false;
      draft.startCoordinate = null;
      draft.orientationLabel = 'unchanged';
      draft.displayedStrandLabel = '';
      return null;
    }
  };

  const open = ({ feature, featureLabel = '' }) => {
    draftGeneration += 1;
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
    draftGeneration += 1;
    Object.assign(draft, initialDraft());
  };

  // Opening Record actions is the explicit record action that reads records a
  // Session Load left unread (776a2f93). It reads once, then recomputes.
  const readRecords = async () => {
    if (!draft.active || draft.reading || draft.pending) return null;
    recompute();
    if (draft.targetFailure !== RECORD_TARGET_NOT_DISCOVERED) return null;
    const generation = draftGeneration;
    Object.assign(draft, {
      reading: true, disabledReason: '', status: 'Reading records…', statusKind: 'pending'
    });
    let outcome;
    try {
      outcome = await action.readRecords();
    } catch (error) {
      outcome = { status: 'error', reason: error?.message || String(error) };
    }
    if (generation !== draftGeneration) return outcome;
    Object.assign(draft, { reading: false, status: '', statusKind: 'idle' });
    recompute();
    if (draft.targetFailure === RECORD_TARGET_NOT_DISCOVERED && outcome?.reason) {
      draft.disabledReason = outcome.reason;
    }
    return outcome;
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
    readRecords,
    setAnchor,
    setOffset,
    setOrientForward,
    placeAtFeatureEnd,
    apply
  });
};
