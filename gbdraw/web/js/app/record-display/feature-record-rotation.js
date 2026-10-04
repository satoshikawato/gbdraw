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

const emptyCapability = () => ({ enabled: false, code: null, message: '' });

const copyCapability = (capability) => ({
  enabled: Boolean(capability?.enabled),
  code: capability?.code || null,
  message: String(capability?.message || '')
});

// "Put this feature at": each choice names an outcome. Custom exposes the raw
// rule "the record starts at a reference point of this feature, shifted by an
// offset"; its references are the existing anchors and the feature-end boundary.
const POSITION_CHOICES = Object.freeze([
  ['start', 'Start of the record'],
  ['end', 'End of the record'],
  ['custom', 'Custom position']
]);
const REFERENCE_CHOICES = Object.freeze([
  ['five-prime', 'this feature’s 5′ end'],
  ['midpoint', 'this feature’s midpoint'],
  ['three-prime', 'this feature’s 3′ end'],
  ['feature-end', 'just after the feature']
]);

const choicesFrom = (entries, capabilityFor = () => emptyCapability()) => entries
  .map(([value, label]) => ({ value, label, ...copyCapability(capabilityFor(value)) }));

const referenceCapability = (capabilities, reference) => (reference === 'feature-end'
  ? capabilities.featureEnd
  : capabilities.anchors?.[reference]);

const positionCapability = (capabilities, position) => {
  if (position === 'start') return capabilities.featureStart;
  if (position === 'end') return capabilities.featureEnd;
  const references = REFERENCE_CHOICES.map(([reference]) => (
    referenceCapability(capabilities, reference)
  ));
  return references.find((capability) => capability?.enabled) || references[0];
};

// A selection stays valid: an unavailable or unknown choice falls back to the
// first available choice in display order.
const validChoice = (choices, selected) => (
  choices.find(({ value }) => value === selected)?.enabled
    ? selected
    : choices.find(({ enabled }) => enabled)?.value ?? selected
);

// The radio choice maps to one resolver intent; the domain math stays in
// resolveFeatureAnchor.
const intentFor = ({ position, reference, offsetText, orientForward }) => {
  const oriented = Boolean(orientForward);
  if (position !== 'custom') {
    return {
      placement: position === 'end' ? 'feature-end' : 'feature-start',
      anchor: null,
      offsetBp: 0,
      orientForward: oriented
    };
  }
  const offsetBp = parseSignedOffset(offsetText).value;
  return reference === 'feature-end'
    ? { placement: 'feature-end', anchor: null, offsetBp, orientForward: oriented }
    : { placement: 'anchor', anchor: reference, offsetBp, orientForward: oriented };
};

const STRAND_SYMBOLS = Object.freeze({ '+': '+', '-': '−' });

const previewFor = ({ recordLabel, resolved, committedReverseComplement }) => {
  if (!Number.isSafeInteger(resolved.startCoordinate)) return '';
  const orientation = resolved.reverseComplement === committedReverseComplement
    ? 'orientation unchanged'
    : resolved.reverseComplement
      ? 'reverse-complemented'
      : 'no longer reverse-complemented';
  const before = STRAND_SYMBOLS[resolved.displayedStrand?.before];
  const after = STRAND_SYMBOLS[resolved.displayedStrand?.after];
  const strandChange = before && after && before !== after
    ? ` · feature strand ${before} → ${after}`
    : '';
  return `${recordLabel} will start at ${resolved.startCoordinate.toLocaleString('en-US')}`
    + ` · ${orientation}${strandChange}`;
};

const initialDraft = () => ({
  active: false,
  feature: null,
  identity: { recordKey: '', biologicalFeatureId: '' },
  recordLabel: '',
  position: 'start',
  reference: 'five-prime',
  offsetText: '0',
  orientForward: false,
  positionChoices: choicesFrom(POSITION_CHOICES),
  referenceChoices: choicesFrom(REFERENCE_CHOICES),
  orientCapability: emptyCapability(),
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
  preview: ''
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

  const resolveDraft = () => action.resolve({ feature: draft.feature, intent: intentFor(draft) });

  const recompute = () => {
    if (!draft.active || !draft.feature) return null;
    try {
      let snapshot = resolveDraft();
      draft.feature = snapshot.feature;
      draft.identity = featureIdentity(snapshot.feature);
      draft.recordLabel = String(snapshot.target.recordId || snapshot.target.recordKey);
      const { capabilities } = snapshot.resolved;
      const placements = [
        capabilities.featureStart,
        ...Object.values(capabilities.anchors || {}),
        capabilities.featureEnd
      ];
      // A limit that disables every placement (for example, a linear record)
      // is a whole-target limit: the reason line states it once (F4).
      const targetLimit = placements.every((capability) => !capability?.enabled)
        ? placements.find((capability) => capability?.message) || placements[0] || emptyCapability()
        : null;
      draft.targetFailure = targetLimit ? String(targetLimit.code || 'unavailable') : '';
      const controlCapability = (capability) => ({
        ...copyCapability(capability),
        ...(targetLimit ? { message: '' } : {})
      });
      draft.positionChoices = choicesFrom(POSITION_CHOICES, (position) => (
        controlCapability(positionCapability(capabilities, position))
      ));
      draft.referenceChoices = choicesFrom(REFERENCE_CHOICES, (reference) => (
        controlCapability(referenceCapability(capabilities, reference))
      ));
      draft.orientCapability = controlCapability(capabilities.orientForward);
      if (!targetLimit) {
        const position = validChoice(draft.positionChoices, draft.position);
        const reference = validChoice(draft.referenceChoices, draft.reference);
        const orientForward = draft.orientForward && draft.orientCapability.enabled;
        if (position !== draft.position || reference !== draft.reference
          || orientForward !== draft.orientForward) {
          Object.assign(draft, { position, reference, orientForward });
          snapshot = resolveDraft();
        }
      }
      // The offset field exists only for Custom; its error is shown next to it.
      const offset = draft.position === 'custom' && !targetLimit
        ? parseSignedOffset(draft.offsetText)
        : { valid: true, message: '' };
      const { eligibility } = snapshot.resolved;
      draft.offsetError = offset.message;
      draft.disabledReason = targetLimit?.message || (offset.valid ? eligibility.message || '' : '');
      draft.canApply = offset.valid && eligibility.enabled && !draft.pending;
      draft.startCoordinate = snapshot.resolved.startCoordinate;
      draft.preview = previewFor({
        recordLabel: draft.recordLabel,
        resolved: snapshot.resolved,
        committedReverseComplement: snapshot.target.committedReverseComplement
      });
      return snapshot;
    } catch (error) {
      // A whole-target failure has one reason; per-control messages are only
      // for control-specific limits of a resolved target.
      draft.targetFailure = error?.kind || 'unavailable';
      draft.positionChoices = choicesFrom(POSITION_CHOICES);
      draft.referenceChoices = choicesFrom(REFERENCE_CHOICES);
      draft.orientCapability = emptyCapability();
      draft.offsetError = '';
      draft.disabledReason = error.message;
      draft.canApply = false;
      draft.startCoordinate = null;
      draft.preview = '';
      return null;
    }
  };

  const open = ({ feature }) => {
    draftGeneration += 1;
    Object.assign(draft, initialDraft(), {
      active: true,
      feature,
      identity: featureIdentity(feature)
    });
    recompute();
    return draft;
  };

  const close = () => {
    draftGeneration += 1;
    Object.assign(draft, initialDraft());
  };

  // Opening the rotation section is the explicit record action that reads
  // records a Session Load left unread (776a2f93). It reads once, then
  // recomputes.
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

  const setPosition = (position) => {
    draft.position = String(position || '');
    recompute();
  };

  const setReference = (reference) => {
    draft.reference = String(reference || '');
    recompute();
  };

  const setOffset = (value) => {
    draft.offsetText = String(value ?? '');
    recompute();
  };

  const setOrientForward = (value) => {
    draft.orientForward = Boolean(value);
    recompute();
  };

  const apply = async () => {
    if (draft.pending) return { status: 'pending' };
    const snapshot = recompute();
    if (!snapshot || !draft.canApply) {
      draft.status = draft.disabledReason || draft.offsetError || 'Record rotation is unavailable.';
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
        intent: intentFor(draft)
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
    setPosition,
    setReference,
    setOffset,
    setOrientForward,
    apply
  });
};
