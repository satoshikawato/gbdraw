// @ts-check
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
 * The popup feature as the Result lists it; its fields come from the rendered
 * feature records (Python owns the shape, R7).
 * @typedef {Record<string, any>} PopupFeature
 */

/**
 * A resolved record transform that Generate or a candidate applies.
 * @typedef {object} ResolvedRecordTransform
 * @property {number} startCoordinate
 * @property {boolean} reverseComplement
 * @property {Record<string, any> | null} anchorIntent
 */

/**
 * @typedef {object} FeatureRecordRotationActionOptions
 * @property {(feature: PopupFeature) => { row: Record<string, any>, target: Record<string, any> }} targetForFeature
 *   record display port: the row and target of the feature's record
 * @property {(row: Record<string, any>, transform: ResolvedRecordTransform) => any} setResolvedTransform
 *   record display port: stages the transform in the record's draft (Apply on Generate)
 * @property {() => Record<string, any> | null} getCommittedSession
 * @property {(args: { committed: Record<string, any> | null, target: Record<string, any>, transform: Record<string, any> }) => { canonical: Record<string, any>, receipt: any }} projectCommittedRecordTransform
 * @property {(run: { canonical: Record<string, any>, row: Record<string, any>, transform: ResolvedRecordTransform }) => Promise<Record<string, any>>} runRecordRotation
 *   composition root port: runs the candidate and commits the target draft
 * @property {((identity: { recordKey: string, biologicalFeatureId: string }) => PopupFeature | null) | null} [resolveCurrentFeature]
 * @property {((feature: PopupFeature) => boolean) | null} [isCurrentFeature]
 * @property {(() => Promise<Record<string, any>>) | null} [readRecords]
 */

/**
 * The action that the workflow drives; `createFeatureRecordRotationAction` returns it.
 * @typedef {object} FeatureRecordRotationActionPort
 * @property {(feature: PopupFeature) => PopupFeature} currentFeature
 * @property {(args: { feature: PopupFeature, intent: Record<string, any> }) => Record<string, any>} resolve
 * @property {(args: { feature: PopupFeature, intent: Record<string, any> }) => Promise<Record<string, any>>} apply
 * @property {(args: { feature: PopupFeature, intent: Record<string, any> }) => Promise<Record<string, any>>} stage
 * @property {() => Promise<Record<string, any>>} readRecords
 */

/**
 * Coordinate one explicit popup feature through the existing domain owners.
 * The controller owns no request editing, generation, admission, or History stack.
 * It receives ports (R13): `targetForFeature` and `setResolvedTransform` of the
 * record display, and `runRecordRotation({ canonical, row, transform })`, which
 * the composition root wires to run the candidate and commit the target draft.
 * @param {FeatureRecordRotationActionOptions} options
 * @returns {Readonly<FeatureRecordRotationActionPort>}
 */
export const createFeatureRecordRotationAction = ({
  targetForFeature,
  setResolvedTransform,
  getCommittedSession,
  projectCommittedRecordTransform,
  runRecordRotation,
  resolveCurrentFeature = null,
  isCurrentFeature = null,
  readRecords = null
}) => {
  if (typeof targetForFeature !== 'function'
    || typeof getCommittedSession !== 'function'
    || typeof projectCommittedRecordTransform !== 'function'
    || typeof runRecordRotation !== 'function') {
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
    // The bound input can differ before the popup opens (a File replaced
    // after Generate, or a Session that needs one Generate), not only after.
    if (typeof isCurrentFeature === 'function' && !isCurrentFeature(resolved)) {
      throw new Error("This feature's input file differs from the one the current Result was drawn from. "
        + 'Generate Diagram to use the current file.');
    }
    return resolved;
  };

  const resolve = ({ feature, intent }) => {
    const resolvedFeature = currentFeature(feature);
    const { row, target } = targetForFeature(resolvedFeature);
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

  const resolveEligible = ({ feature, intent }) => {
    const snapshot = resolve({ feature, intent });
    if (!snapshot.resolved.eligibility.enabled) {
      throw new Error(snapshot.resolved.eligibility.message);
    }
    return snapshot;
  };

  const resolvedTransform = (resolved) => ({
    startCoordinate: resolved.startCoordinate,
    reverseComplement: resolved.reverseComplement,
    anchorIntent: resolved.provenance
  });

  const apply = async ({ feature, intent }) => {
    const { row, target, resolved } = resolveEligible({ feature, intent });
    const projection = projectCommittedRecordTransform({
      committed: getCommittedSession(),
      target,
      transform: {
        recordLength: target.recordLength,
        startCoordinate: resolved.startCoordinate,
        reverseComplement: resolved.reverseComplement
      }
    });
    const outcome = await runRecordRotation({
      canonical: projection.canonical,
      row,
      transform: resolvedTransform(resolved)
    });
    return { ...outcome, receipt: projection.receipt, resolved };
  };

  // Apply on Generate (PD-OI-085): the same resolved transform goes into the
  // record display draft that the sidebar edits and Generate applies, as one
  // undoable step. No candidate runs and the Result is untouched.
  const stage = async ({ feature, intent }) => {
    const { row, resolved } = resolveEligible({ feature, intent });
    const busy = await setResolvedTransform(row, resolvedTransform(resolved));
    return busy?.status === 'busy' ? busy : { status: 'ok', resolved };
  };

  // The explicit record action that reads records a Session Load left unread.
  const readTargetRecords = async () => (typeof readRecords === 'function'
    ? readRecords()
    : { status: 'unavailable', reason: 'Records cannot be read here.' });

  return Object.freeze({ currentFeature, resolve, apply, stage, readRecords: readTargetRecords });
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

/**
 * @param {readonly string[][]} entries `[value, label]` pairs
 * @param {(value: string) => Record<string, any> | undefined} [capabilityFor]
 */
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

// "<record> will start at <n> · <orientation>", relative to the Result.
const startSentence = ({ recordLabel, startCoordinate, reverseComplement, committedReverseComplement }) => {
  const orientation = reverseComplement === committedReverseComplement
    ? 'orientation unchanged'
    : reverseComplement
      ? 'reverse-complemented'
      : 'no longer reverse-complemented';
  return `${recordLabel} will start at ${startCoordinate.toLocaleString('en-US')} · ${orientation}`;
};

const previewFor = ({ recordLabel, resolved, committedReverseComplement }) => {
  if (!Number.isSafeInteger(resolved.startCoordinate)) return '';
  const before = STRAND_SYMBOLS[resolved.displayedStrand?.before];
  const after = STRAND_SYMBOLS[resolved.displayedStrand?.after];
  const strandChange = before && after && before !== after
    ? ` · feature strand ${before} → ${after}`
    : '';
  return startSentence({
    recordLabel,
    startCoordinate: resolved.startCoordinate,
    reverseComplement: resolved.reverseComplement,
    committedReverseComplement
  }) + strandChange;
};

// The start and orientation staged in the record's draft for Generate, when
// they differ from the Result. A draft without a start keeps the default start:
// the record's first base, or its last base when reverse-complemented.
const pendingStartFor = (target) => {
  const pending = target.pendingTransform;
  if (!pending) return null;
  const startCoordinate = pending.startCoordinate
    ?? (pending.reverseComplement ? target.recordLength : 1);
  return Number.isSafeInteger(startCoordinate)
    ? { startCoordinate, reverseComplement: pending.reverseComplement }
    : null;
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
  preview: '',
  // The record's start and orientation already staged for Generate (PD-OI-085).
  pendingPreview: ''
});

/**
 * @typedef {object} FeatureRecordRotationWorkflowOptions
 * @property {FeatureRecordRotationActionPort} action
 * @property {(<T extends object>(value: T) => T) | undefined} [makeReactive] Vue `reactive`
 * @property {((feature: PopupFeature, identity: { recordKey: string, biologicalFeatureId: string }) => void) | null} [onRebind]
 */

/**
 * Own only the popup's ephemeral draft and presentation state. Domain math and
 * transactional generation stay in createFeatureRecordRotationAction.
 * @param {FeatureRecordRotationWorkflowOptions} options
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
      const { committedReverseComplement } = snapshot.target;
      const pending = pendingStartFor(snapshot.target);
      draft.pendingPreview = pending
        ? `Pending for Generate: ${startSentence({
          recordLabel: draft.recordLabel, ...pending, committedReverseComplement
        })}`
        : '';
      // A choice that resolves to the staged transform is stated once, by the pending line.
      draft.preview = pending
        && pending.startCoordinate === snapshot.resolved.startCoordinate
        && pending.reverseComplement === snapshot.resolved.reverseComplement
        ? ''
        : previewFor({
          recordLabel: draft.recordLabel,
          resolved: snapshot.resolved,
          committedReverseComplement
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
      draft.pendingPreview = '';
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

  // Both apply buttons share this validation (canApply). A failure keeps its
  // one reason, on the reason line or next to the offset field (F4).
  const applicableSnapshot = () => {
    const snapshot = recompute();
    if (snapshot && draft.canApply) return snapshot;
    const shown = Boolean(draft.disabledReason || draft.offsetError);
    draft.status = shown ? '' : 'Record rotation is unavailable.';
    draft.statusKind = shown ? 'idle' : 'error';
    return null;
  };

  // A failure message for the status line, unless the reason line, recomputed
  // after the action, already states it.
  const unlessShown = (message, fallback) => (
    message && message !== draft.disabledReason ? message : fallback
  );

  const apply = async () => {
    if (draft.pending) return { status: 'pending' };
    const snapshot = applicableSnapshot();
    if (!snapshot) return { status: 'disabled' };
    // F7: a popup closed or retargeted meanwhile has a new draft; this run
    // writes nothing into it.
    const generation = draftGeneration;
    draft.pending = true;
    draft.canApply = false;
    draft.status = 'Regenerating the target record…';
    draft.statusKind = 'pending';
    let outcome;
    let failure = '';
    try {
      outcome = await action.apply({
        feature: snapshot.feature,
        intent: intentFor(draft)
      });
      if (generation === draftGeneration && outcome?.status === 'ok') {
        const rebound = action.currentFeature(snapshot.feature);
        draft.feature = rebound;
        if (typeof onRebind === 'function') onRebind(rebound, draft.identity);
      }
    } catch (error) {
      outcome = { status: 'error', error };
      failure = error?.message || String(error);
    }
    if (generation !== draftGeneration) return outcome;
    draft.pending = false;
    recompute();
    const kept = 'The previous Result was kept.';
    if (outcome?.status === 'ok') {
      draft.status = 'Record rotation applied and regenerated.';
      draft.statusKind = 'success';
    } else if (outcome?.status === 'canceled') {
      draft.status = `Record rotation canceled. ${kept}`;
      draft.statusKind = 'idle';
    } else {
      draft.status = `${unlessShown(failure, 'Record rotation failed.')} ${kept}`;
      draft.statusKind = 'error';
    }
    return outcome;
  };

  // Apply on Generate: stage the resolved transform; Generate Diagram applies
  // every staged record together (PD-OI-085).
  let staging = false;
  const stage = async () => {
    if (draft.pending || staging) return { status: 'pending' };
    const snapshot = applicableSnapshot();
    if (!snapshot) return { status: 'disabled' };
    const generation = draftGeneration;
    staging = true;
    let outcome;
    try {
      outcome = await action.stage({
        feature: snapshot.feature,
        intent: intentFor(draft)
      });
    } catch (error) {
      outcome = { status: 'error', error };
    } finally {
      staging = false;
    }
    if (generation !== draftGeneration) return outcome;
    recompute();
    if (outcome?.status === 'ok') {
      draft.status = 'Record rotation will apply on the next Generate Diagram.';
      draft.statusKind = 'success';
    } else {
      draft.status = outcome?.status === 'busy'
        ? outcome.reason
        : unlessShown(outcome?.error?.message, '');
      draft.statusKind = draft.status ? 'error' : 'idle';
    }
    return outcome;
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
    apply,
    stage
  });
};
