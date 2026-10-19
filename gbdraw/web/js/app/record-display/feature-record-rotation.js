// @ts-check
import { resolveFeatureAnchor } from './feature-anchor.js';
import { RECORD_TARGET_NOT_DISCOVERED } from '../record-display-options.js';
import { msg, t, tm } from '../../utils/i18n.js';
/** @import { CommittedDisplayRow, ResolvedRecordTransform } from '../record-display-options.js' */
/** @import { AnchorIntent } from './feature-anchor.js' */
/** @import { SessionBusy } from '../../state.js' */
/** @import { RecordKeyedFeature } from '../../services/feature-identity.js' */
/** @import { CommittedCanonicalSession, RecordTransformInput, RecordTransformTarget } from '../../services/session-request.js' */

/**
 * One capability verdict as the popup shows it.
 * @typedef {{ enabled: boolean, code: string | null, message: string }} CapabilityView
 */

/** @param {Record<string, any> | null | undefined} feature */
const featureIdentity = (feature) => ({
  recordKey: String(feature?.record_key ?? feature?.recordKey ?? ''),
  biologicalFeatureId: String(
    feature?.biological_feature_id ?? feature?.biologicalFeatureId ?? ''
  )
});

/** @param {Record<string, any> | null | undefined} feature @param {{ recordKey: string, biologicalFeatureId: string }} identity */
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
 * @typedef {object} FeatureRecordRotationActionOptions
 * @property {(feature: PopupFeature) => { row: CommittedDisplayRow, target: Record<string, any> }} targetForFeature
 *   record display port: the row and target of the feature's record
 * @property {(row: CommittedDisplayRow, transform: ResolvedRecordTransform) => any} setResolvedTransform
 *   record display port: stages the transform in the record's draft (Apply on Generate)
 * @property {() => CommittedCanonicalSession | null} getCommittedSession
 * @property {(args: { committed: CommittedCanonicalSession | null, target: RecordTransformTarget, transform: RecordTransformInput }) => { canonical: Record<string, any>, receipt: any }} projectCommittedRecordTransform
 * @property {(run: { canonical: Record<string, any>, row: CommittedDisplayRow, transform: ResolvedRecordTransform }) => Promise<Record<string, any>>} runRecordRotation
 *   composition root port: runs the candidate and commits the target draft
 * @property {((identity: { recordKey: string, biologicalFeatureId: string }) => PopupFeature | null) | null} [resolveCurrentFeature]
 * @property {((feature: RecordKeyedFeature) => boolean) | null} [isCurrentFeature]
 * @property {(() => Promise<Record<string, any>> | SessionBusy) | null} [readRecords]
 */

/**
 * The action that the workflow drives; `createFeatureRecordRotationAction` returns it.
 * @typedef {object} FeatureRecordRotationActionPort
 * @property {(feature: PopupFeature) => PopupFeature} currentFeature
 * @property {(args: { feature: PopupFeature, intent: AnchorIntent }) => Record<string, any>} resolve
 * @property {(args: { feature: PopupFeature, intent: AnchorIntent }) => Promise<Record<string, any>>} apply
 * @property {(args: { feature: PopupFeature, intent: AnchorIntent }) => Promise<Record<string, any>>} stage
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

  /** @param {PopupFeature} feature @returns {PopupFeature} */
  const currentFeature = (feature) => {
    const identity = featureIdentity(feature);
    if (!identity.recordKey || !identity.biologicalFeatureId) {
      throw new Error('Feature rotation requires a stable record and feature identity.');
    }
    const resolved = typeof resolveCurrentFeature === 'function'
      ? resolveCurrentFeature(identity)
      : feature;
    if (!resolved || !sameFeatureIdentity(resolved, identity)) {
      throw new Error(t('The popup feature is no longer present in the current Result.'));
    }
    // The bound input can differ before the popup opens (a File replaced
    // after Generate, or a Session that needs one Generate), not only after.
    if (typeof isCurrentFeature === 'function' && !isCurrentFeature(resolved)) {
      throw new Error(t("This feature's input file differs from the one the current Result was drawn from.")
        + ' ' + t('Generate Diagram to use the current file.'));
    }
    return resolved;
  };

  /** @param {{ feature: PopupFeature, intent: AnchorIntent }} args */
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

  /** @param {{ feature: PopupFeature, intent: AnchorIntent }} args */
  const resolveEligible = ({ feature, intent }) => {
    const snapshot = resolve({ feature, intent });
    if (!snapshot.resolved.eligibility.enabled) {
      throw new Error(snapshot.resolved.eligibility.message);
    }
    return snapshot;
  };

  // `resolveEligible` throws unless the resolution is eligible, and an
  // eligible resolution sets the start coordinate, the orientation and the
  // provenance (`resolveFeatureAnchor`, feature-anchor.js).
  /** @param {{ startCoordinate: number | null, reverseComplement: boolean | null, provenance: Record<string, unknown> | null }} resolved */
  const resolvedTransform = (resolved) => /** @type {ResolvedRecordTransform} */ ({
    startCoordinate: resolved.startCoordinate,
    reverseComplement: resolved.reverseComplement,
    anchorIntent: resolved.provenance
  });

  /** @param {{ feature: PopupFeature, intent: AnchorIntent }} args */
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
  /** @param {{ feature: PopupFeature, intent: AnchorIntent }} args */
  const stage = async ({ feature, intent }) => {
    const { row, resolved } = resolveEligible({ feature, intent });
    const busy = await setResolvedTransform(row, resolvedTransform(resolved));
    return busy?.status === 'busy' ? busy : { status: 'ok', resolved };
  };

  // The explicit record action that reads records a Session Load left unread.
  const readTargetRecords = async () => (typeof readRecords === 'function'
    ? readRecords()
    : { status: 'unavailable', reason: t('Records cannot be read here.') });

  return Object.freeze({ currentFeature, resolve, apply, stage, readRecords: readTargetRecords });
};

/** @param {unknown} value */
const parseSignedOffset = (value) => {
  const text = String(value ?? '').trim();
  if (!text) {
    return { valid: false, value: Number.NaN, message: t('Enter a signed whole-number offset.') };
  }
  if (!/^[+-]?\d+$/.test(text)) {
    return {
      valid: false,
      value: Number.NaN,
      message: t('Offset must be a signed whole number of base pairs.')
    };
  }
  const number = Number(text);
  if (!Number.isSafeInteger(number)) {
    return {
      valid: false,
      value: Number.NaN,
      message: t('Offset is outside the supported integer range.')
    };
  }
  return { valid: true, value: number, message: '' };
};

/** @returns {CapabilityView} */
const emptyCapability = () => ({ enabled: false, code: null, message: '' });

/**
 * @param {{ enabled?: unknown, code?: string | null, message?: unknown } | null | undefined} capability
 * @returns {CapabilityView}
 */
const copyCapability = (capability) => ({
  enabled: Boolean(capability?.enabled),
  code: capability?.code || null,
  message: String(capability?.message || '')
});

// "Put this feature at": each choice names an outcome. Custom exposes the raw
// rule "the record starts at a reference point of this feature, shifted by an
// offset"; its references are the existing anchors and the feature-end boundary.
const POSITION_CHOICES = Object.freeze([
  ['start', msg('Start of the record')],
  ['end', msg('End of the record')],
  ['custom', msg('Custom position')]
]);
const REFERENCE_CHOICES = Object.freeze([
  ['five-prime', msg('this feature’s 5′ end')],
  ['midpoint', msg('this feature’s midpoint')],
  ['three-prime', msg('this feature’s 3′ end')],
  ['feature-end', msg('just after the feature')]
]);

/**
 * @param {readonly string[][]} entries `[value, label]` pairs
 * @param {(value: string) => CapabilityView | undefined} [capabilityFor]
 */
const choicesFrom = (entries, capabilityFor = () => emptyCapability()) => entries
  .map(([value, label]) => ({ value, label: tm(label), ...copyCapability(capabilityFor(value)) }));

/**
 * @param {{ featureEnd: CapabilityView, anchors?: Readonly<Record<string, CapabilityView | undefined>> }} capabilities
 * @param {string} reference One of REFERENCE_CHOICES.
 */
const referenceCapability = (capabilities, reference) => (reference === 'feature-end'
  ? capabilities.featureEnd
  : capabilities.anchors?.[reference]);

/**
 * @param {{ featureStart: CapabilityView, featureEnd: CapabilityView, anchors?: Readonly<Record<string, CapabilityView | undefined>> }} capabilities
 * @param {string} position One of POSITION_CHOICES.
 */
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
/** @param {{ value: string, enabled: boolean }[]} choices @param {string} selected */
const validChoice = (choices, selected) => (
  choices.find(({ value }) => value === selected)?.enabled
    ? selected
    : choices.find(({ enabled }) => enabled)?.value ?? selected
);

// The radio choice maps to one resolver intent; the domain math stays in
// resolveFeatureAnchor.
/**
 * @param {{ position: string, reference: string, offsetText: string, orientForward: boolean }} draft
 * @returns {AnchorIntent}
 */
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

// The strand is '+', '-', or another value of the resolver's validated strand
// domain ('unstranded', 'mixed'), or null; none names an Object.prototype member.
/** @type {Readonly<Record<string, string | undefined>>} */
const STRAND_SYMBOLS = Object.freeze({ '+': '+', '-': '−' });

// "<record> will start at <n> · <orientation>", relative to the Result.
/** @param {{ recordLabel: string, startCoordinate: number, reverseComplement: boolean | null, committedReverseComplement: boolean }} sentence */
const startSentence = ({ recordLabel, startCoordinate, reverseComplement, committedReverseComplement }) => {
  const orientation = reverseComplement === committedReverseComplement
    ? t('orientation unchanged')
    : reverseComplement
      ? t('reverse-complemented')
      : t('no longer reverse-complemented');
  return `${t('{record} will start at {start}', { record: recordLabel, start: startCoordinate.toLocaleString('en-US') })} · ${orientation}`;
};

/**
 * @param {{
 *   recordLabel: string, committedReverseComplement: boolean,
 *   resolved: { startCoordinate: number | null, reverseComplement: boolean | null, displayedStrand?: { before: string | null, after: string | null } }
 * }} preview
 */
const previewFor = ({ recordLabel, resolved, committedReverseComplement }) => {
  if (!Number.isSafeInteger(resolved.startCoordinate)) return '';
  // A safe-integer start means an eligible resolution, which sets both strands.
  const before = STRAND_SYMBOLS[/** @type {string} */ (resolved.displayedStrand?.before)];
  const after = STRAND_SYMBOLS[/** @type {string} */ (resolved.displayedStrand?.after)];
  const strandChange = before && after && before !== after
    ? ` · ${t('feature strand {before} → {after}', { before, after })}`
    : '';
  return startSentence({
    recordLabel,
    // `isSafeInteger` above is true only for a number.
    startCoordinate: /** @type {number} */ (resolved.startCoordinate),
    reverseComplement: resolved.reverseComplement,
    committedReverseComplement
  }) + strandChange;
};

// The start and orientation staged in the record's draft for Generate, when
// they differ from the Result. A draft without a start keeps the default start:
// the record's first base, or its last base when reverse-complemented.
/** @param {{ recordLength: number | null | undefined, pendingTransform?: { startCoordinate: number | null, reverseComplement: boolean } | null }} target */
const pendingStartFor = (target) => {
  const pending = target.pendingTransform;
  if (!pending) return null;
  const startCoordinate = pending.startCoordinate
    ?? (pending.reverseComplement ? target.recordLength : 1);
  return Number.isSafeInteger(startCoordinate)
    // `isSafeInteger` is true only for a number.
    ? { startCoordinate: /** @type {number} */ (startCoordinate), reverseComplement: pending.reverseComplement }
    : null;
};

const initialDraft = () => ({
  active: false,
  feature: /** @type {PopupFeature | null} */ (null),
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
  startCoordinate: /** @type {number | null} */ (null),
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

  // recompute is the only caller; it returns early without a feature and stores the snapshot's feature before its second call.
  const resolveDraft = () => action.resolve({ feature: /** @type {PopupFeature} */ (draft.feature), intent: intentFor(draft) });

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
      const targetLimit = placements.every((/** @type {CapabilityView | undefined} */ capability) => !capability?.enabled)
        ? placements.find((/** @type {CapabilityView | undefined} */ capability) => capability?.message) || placements[0] || emptyCapability()
        : null;
      draft.targetFailure = targetLimit ? String(targetLimit.code || 'unavailable') : '';
      /** @param {{ enabled?: unknown, code?: string | null, message?: unknown } | null | undefined} capability */
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
        ? t('Pending for Generate: {sentence}', { sentence: startSentence({
          recordLabel: draft.recordLabel, ...pending, committedReverseComplement
        }) })
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

  /** @param {{ feature: PopupFeature }} args */
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
      reading: true, disabledReason: '', status: t('Reading records…'), statusKind: 'pending'
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

  /** @param {unknown} position */
  const setPosition = (position) => {
    draft.position = String(position || '');
    recompute();
  };

  /** @param {unknown} reference */
  const setReference = (reference) => {
    draft.reference = String(reference || '');
    recompute();
  };

  /** @param {unknown} value */
  const setOffset = (value) => {
    draft.offsetText = String(value ?? '');
    recompute();
  };

  /** @param {unknown} value */
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
    draft.status = shown ? '' : t('Record rotation is unavailable.');
    draft.statusKind = shown ? 'idle' : 'error';
    return null;
  };

  // A failure message for the status line, unless the reason line, recomputed
  // after the action, already states it.
  /** @param {string | undefined} message @param {string} fallback */
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
    draft.status = t('Regenerating the target record…');
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
    const kept = t('The previous Result was kept.');
    if (outcome?.status === 'ok') {
      draft.status = t('Record rotation applied and regenerated.');
      draft.statusKind = 'success';
    } else if (outcome?.status === 'canceled') {
      draft.status = `${t('Record rotation canceled.')} ${kept}`;
      draft.statusKind = 'idle';
    } else {
      draft.status = `${unlessShown(failure, t('Record rotation failed.'))} ${kept}`;
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
      draft.status = t('Record rotation will apply on the next Generate Diagram.');
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
