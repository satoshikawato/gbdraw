// @ts-check
import { locationPartBound } from '../../services/feature-utils.js';

/**
 * One capability verdict of the anchor resolver. `basis` is set on the offset
 * verdict only.
 * @typedef {object} AnchorCapability
 * @property {boolean} enabled
 * @property {string | null} code
 * @property {string} message
 * @property {string} [basis]
 */

/**
 * One source-location part of a feature, as stored in the extracted feature
 * row. The resolver reads and validates every field, so the values are not
 * trusted.
 * @typedef {object} AnchorPartInput
 * @property {unknown} start
 * @property {unknown} end
 * @property {string} strand
 */

/**
 * A part as the resolver read it: `start` and `end` come from the shared
 * location-part reader (`locationPartBound`, NaN when missing), so the anchor
 * reads a part by the same rule as its anchor profile (OV-1016).
 * @typedef {object} AnchorPart
 * @property {number} start
 * @property {number} end
 * @property {string} strand
 */

/**
 * The source-location profile of a feature. The resolver validates every
 * field against its allowed values, so the fields are plain strings.
 * @typedef {object} AnchorProfile
 * @property {string} precision
 * @property {string} operator
 * @property {string} partOrder
 * @property {string} strand
 */

/**
 * @typedef {object} AnchorIntent
 * @property {string} placement
 * @property {string | null} anchor
 * @property {number} offsetBp
 * @property {boolean} orientForward
 */

/**
 * @typedef {object} AnchorIdentity
 * @property {string} recordKey
 * @property {string} biologicalFeatureId
 */

/**
 * @typedef {object} AnchorRequest
 * @property {number} recordLength
 * @property {boolean} effectiveCircular
 * @property {boolean} [cropped]
 * @property {boolean} currentReverseComplement
 * @property {AnchorIdentity} identity
 * @property {AnchorPartInput[] | undefined} parts
 * @property {AnchorProfile | null | undefined} profile
 * @property {AnchorIntent} intent
 */

/**
 * @param {{ basis?: string }} [extra]
 * @returns {AnchorCapability}
 */
const available = (extra = {}) => ({ enabled: true, code: null, message: '', ...extra });
/**
 * @param {string} code
 * @param {string} message
 * @param {{ basis?: string }} [extra]
 * @returns {AnchorCapability}
 */
const unavailable = (code, message, extra = {}) => ({ enabled: false, code, message, ...extra });

/**
 * @param {AnchorCapability[]} operations
 * @returns {AnchorCapability | null}
 */
const firstUnavailable = (...operations) => operations.find((operation) => !operation.enabled) || null;

/**
 * @param {number} coordinate
 * @param {number} length
 * @returns {number}
 */
const wrapCoordinate = (coordinate, length) => 1 + ((((coordinate - 1) % length) + length) % length);

/**
 * @param {string} sourceStrand
 * @param {boolean} reverseComplement
 * @returns {string}
 */
const displayedStrand = (sourceStrand, reverseComplement) => {
  if (!['+', '-'].includes(sourceStrand)) return sourceStrand;
  if (!reverseComplement) return sourceStrand;
  return sourceStrand === '+' ? '-' : '+';
};

/**
 * @param {Pick<AnchorRequest, 'recordLength' | 'effectiveCircular' | 'cropped' | 'currentReverseComplement'>} options
 * @returns {AnchorCapability}
 */
const recordCapability = ({ recordLength, effectiveCircular, cropped, currentReverseComplement }) => {
  if (!Number.isSafeInteger(recordLength) || recordLength <= 0) {
    return unavailable('record-length-unavailable', 'Record length must be a positive integer.');
  }
  if (effectiveCircular !== true) {
    return unavailable('record-not-circular', 'Record rotation requires an effectively circular record.');
  }
  if (cropped === true) {
    return unavailable('record-cropped', 'Record rotation is unavailable for a cropped record.');
  }
  if (currentReverseComplement !== true && currentReverseComplement !== false) {
    return unavailable('orientation-state-invalid', 'Current record orientation is unavailable.');
  }
  return available();
};

/**
 * @param {AnchorProfile | null | undefined} profile
 * @returns {AnchorCapability}
 */
const profileCapability = (profile) => {
  if (!profile || !['exact', 'fuzzy', 'unavailable'].includes(profile.precision)
    || !['single', 'join', 'order', 'unknown'].includes(profile.operator)
    || !['biological', 'source-forward', 'ambiguous'].includes(profile.partOrder)
    || !['+', '-', 'unstranded', 'mixed'].includes(profile.strand)) {
    return unavailable('location-profile-invalid', 'Feature source-location capabilities are unavailable.');
  }
  if (profile.precision === 'unavailable') {
    return unavailable(
      'feature-metadata-refresh-required',
      'Generate again to refresh feature location metadata before rotating this record.'
    );
  }
  return available();
};

/**
 * @param {string} partStrand
 * @param {string} profileStrand
 * @returns {boolean}
 */
const partStrandMatches = (partStrand, profileStrand) => {
  if (profileStrand === 'mixed') return true;
  if (profileStrand === 'unstranded') {
    return !['+', '-'].includes(partStrand);
  }
  return partStrand === profileStrand;
};

/**
 * Reads each part once; every later check and coordinate uses the read values.
 * @param {AnchorPartInput[] | undefined} parts
 * @returns {Array<AnchorPart | null> | undefined}
 */
const readAnchorParts = (parts) => (Array.isArray(parts)
  ? parts.map((part) => (part ? {
    start: locationPartBound(part, 'start'),
    end: locationPartBound(part, 'end'),
    strand: part.strand
  } : null))
  : undefined);

/**
 * @param {Array<AnchorPart | null> | undefined} parts
 * @param {AnchorProfile} profile
 * @param {number} recordLength
 * @returns {AnchorCapability}
 */
const partsCapability = (parts, profile, recordLength) => {
  if (!Array.isArray(parts) || parts.length === 0) {
    return unavailable('location-parts-unavailable', 'Feature source-location parts are unavailable.');
  }
  for (const part of parts) {
    if (!part || !Number.isSafeInteger(part.start) || !Number.isSafeInteger(part.end)
      || part.start < 0 || part.end <= part.start || part.end > recordLength) {
      return unavailable('location-parts-invalid', 'Feature source-location parts are outside the record.');
    }
    if (!partStrandMatches(part.strand, profile.strand)) {
      return unavailable('location-strand-mismatch', 'Feature part strands disagree with the source profile.');
    }
  }
  // The loop above returned for every missing part.
  const readParts = /** @type {AnchorPart[]} */ (parts);
  const sourceOrder = [...readParts].sort((left, right) => left.start - right.start || left.end - right.end);
  if (sourceOrder.some((part, index) => index > 0 && part.start < sourceOrder[index - 1].end)) {
    return unavailable('location-parts-overlap', 'Overlapping feature parts cannot be traversed safely.');
  }
  if (profile.partOrder === 'source-forward') {
    const wraps = readParts.reduce((count, part, index) => count
      + (index > 0 && part.start < readParts[index - 1].start ? 1 : 0), 0);
    if (wraps > 1) {
      return unavailable('location-source-order-invalid', 'Feature parts do not form one source-forward circular path.');
    }
  }
  return available();
};

/**
 * @param {AnchorCapability} record
 * @param {AnchorCapability} profileState
 * @param {AnchorCapability} partsState
 * @param {AnchorProfile | Record<string, never>} profile - `{}` when no profile exists; the
 *   states above are then unavailable and return first.
 * @returns {AnchorCapability}
 */
const positionalCapability = (record, profileState, partsState, profile) => {
  const common = firstUnavailable(record, profileState, partsState);
  if (common) return common;
  if (profile.precision !== 'exact') {
    return unavailable('location-fuzzy', 'This operation requires exact feature coordinates.');
  }
  if (profile.operator === 'order') {
    return unavailable('location-order-ambiguous', 'Ordered feature parts do not define an exact traversal path.');
  }
  if (profile.operator === 'unknown') {
    return unavailable('location-operator-unknown', 'Feature part traversal semantics are unavailable.');
  }
  if (profile.partOrder === 'ambiguous' || profile.strand === 'mixed') {
    return unavailable('location-traversal-ambiguous', 'Feature part order or strand is ambiguous.');
  }
  return available();
};

/**
 * @param {AnchorPart[]} parts
 * @param {number} offset
 * @param {number} direction
 * @returns {number | null}
 */
const baseAtCoveredOffset = (parts, offset, direction) => {
  let remaining = offset;
  for (const part of parts) {
    const length = part.end - part.start;
    if (remaining < length) {
      return direction === -1 ? part.end - remaining : part.start + remaining + 1;
    }
    remaining -= length;
  }
  return null;
};

/**
 * @param {AnchorIntent} intent
 * @returns {AnchorCapability}
 */
const intentCapability = (intent) => {
  if (!intent || !['feature-start', 'anchor', 'feature-end'].includes(intent.placement)) {
    return unavailable('placement-invalid', 'Choose a feature-start, anchor, or feature-end placement.');
  }
  if (intent.placement === 'anchor'
    && !(/** @type {(string | null)[]} */ (['five-prime', 'midpoint', 'three-prime'])).includes(intent.anchor)) {
    return unavailable('anchor-invalid', 'Choose a supported feature anchor.');
  }
  if (intent.placement !== 'anchor' && intent.anchor !== null) {
    return unavailable('anchor-not-applicable', 'Feature-start and feature-end placements do not use an anchor.');
  }
  if (!Number.isSafeInteger(intent.offsetBp)) {
    return unavailable('offset-invalid', 'Offset must be a safe integer number of base pairs.');
  }
  if (intent.placement === 'feature-start' && intent.offsetBp !== 0) {
    return unavailable('offset-not-applicable', 'Feature-start placement does not use an offset.');
  }
  if (intent.orientForward !== true && intent.orientForward !== false) {
    return unavailable('orient-forward-invalid', 'Orient-forward intent must be true or false.');
  }
  return available();
};

/**
 * @param {AnchorIdentity} identity
 * @returns {AnchorCapability}
 */
const identityCapability = (identity) => {
  if (!identity || typeof identity.recordKey !== 'string' || !identity.recordKey
    || identity.recordKey.includes('\0') || typeof identity.biologicalFeatureId !== 'string'
    || !identity.biologicalFeatureId || identity.biologicalFeatureId.includes('\0')) {
    return unavailable('feature-identity-unavailable', 'Feature rotation requires a stable record and feature identity.');
  }
  return available();
};

/**
 * @param {AnchorRequest} request
 */
export const resolveFeatureAnchor = ({
  recordLength,
  effectiveCircular,
  cropped = false,
  currentReverseComplement,
  identity,
  parts,
  profile,
  intent
}) => {
  const record = recordCapability({ recordLength, effectiveCircular, cropped, currentReverseComplement });
  const profileState = profileCapability(profile);
  const readParts = readAnchorParts(parts);
  const partsState = profileState.enabled
    ? partsCapability(readParts, /** @type {AnchorProfile} */ (profile), recordLength) : profileState;
  const positional = positionalCapability(record, profileState, partsState, profile || {});
  // Each use below is reached only when `positional` (or the eligibility
  // check) found the profile and the parts valid, so both exist there.
  const located = /** @type {AnchorProfile} */ (profile);
  const locatedParts = /** @type {AnchorPart[]} */ (readParts);

  const midpoint = positional;
  const biologicalEnds = !positional.enabled ? positional
    : !['+', '-'].includes(located.strand)
      ? unavailable('biological-direction-unavailable', 'Biological 5-prime and 3-prime ends require a known strand.')
      : available();
  const offsetBasis = profile?.strand === 'unstranded' ? 'source-forward' : 'strand';
  const offset = !positional.enabled ? { ...positional, basis: offsetBasis }
    : available({ basis: offsetBasis });
  const orientForward = !record.enabled ? record
    : !profileState.enabled ? profileState
      : !['+', '-'].includes(located.strand)
        ? unavailable('orientation-strand-unavailable', 'Orient forward requires one known feature strand.')
        : available();
  const featureEnd = !positional.enabled ? positional
    : located.strand === 'unstranded' && located.partOrder !== 'source-forward'
      ? unavailable('feature-end-direction-unavailable', 'Feature-end placement requires a known display traversal path.')
      : available();
  // The feature's first base in display order is its 5-prime or 3-prime base,
  // so placing it at the start needs a known strand.
  const featureStart = !positional.enabled ? positional
    : !['+', '-'].includes(located.strand)
      ? unavailable('feature-start-direction-unavailable', 'Needs a known feature strand.')
      : available();
  const capabilities = {
    featureStart,
    anchors: {
      'five-prime': biologicalEnds,
      midpoint,
      'three-prime': biologicalEnds
    },
    offset,
    orientForward,
    featureEnd
  };

  const identityState = identityCapability(identity);
  const intentState = intentCapability(intent);
  let eligibility = firstUnavailable(identityState, intentState);
  if (!eligibility && intent.orientForward) eligibility = firstUnavailable(orientForward);
  if (!eligibility) eligibility = intent.placement === 'feature-end'
    ? firstUnavailable(featureEnd, offset)
    : intent.placement === 'feature-start'
      ? firstUnavailable(featureStart)
      : firstUnavailable(capabilities.anchors[/** @type {'five-prime' | 'midpoint' | 'three-prime'} */ (intent.anchor)], offset);
  if (eligibility) {
    return {
      eligibility,
      capabilities,
      startCoordinate: null,
      reverseComplement: null,
      sourceAnchorCoordinate: null,
      sourceOutgoingBoundary: null,
      displayedStrand: { before: null, after: null },
      provenance: null
    };
  }

  const reverseComplement = intent.orientForward ? located.strand === '-' : currentReverseComplement;
  const strandAfter = displayedStrand(located.strand, reverseComplement);
  // Feature-start is the anchor whose base comes first in the resulting display
  // order: the 5-prime base when the feature reads forward there, else the
  // 3-prime base. It is recorded as that anchor, so provenance keeps its format.
  const placement = intent.placement === 'feature-start' ? 'anchor' : intent.placement;
  const anchor = intent.placement === 'feature-start'
    ? (strandAfter === '+' ? 'five-prime' : 'three-prime')
    : intent.anchor;
  const direction = located.strand === '-' ? -1 : 1;
  const coveredLength = locatedParts.reduce((total, part) => total + part.end - part.start, 0);
  /** @type {number | null} */
  let sourceAnchorCoordinate = null;
  /** @type {number | null} */
  let sourceOutgoingBoundary = null;
  // Both branches find a base: the parts are valid and the offsets lie inside them.
  /** @type {number} */
  let baseCoordinate;
  if (placement === 'anchor') {
    const coveredOffset = anchor === 'five-prime' ? 0
      : anchor === 'midpoint' ? Math.floor((coveredLength - 1) / 2)
        : coveredLength - 1;
    sourceAnchorCoordinate = baseAtCoveredOffset(locatedParts, coveredOffset, direction);
    baseCoordinate = /** @type {number} */ (sourceAnchorCoordinate);
  } else {
    const displayStep = reverseComplement ? -1 : 1;
    const followsPath = displayStep === direction;
    const lastOffset = followsPath ? coveredLength - 1 : 0;
    const lastIncludedBase = baseAtCoveredOffset(locatedParts, lastOffset, direction);
    sourceOutgoingBoundary = wrapCoordinate(/** @type {number} */ (lastIncludedBase) + displayStep, recordLength);
    baseCoordinate = sourceOutgoingBoundary;
  }
  const offsetDirection = located.strand === '-' ? -1 : 1;
  const startCoordinate = wrapCoordinate(
    baseCoordinate + (offsetDirection * intent.offsetBp),
    recordLength
  );

  return {
    eligibility: available(),
    capabilities,
    startCoordinate,
    reverseComplement,
    sourceAnchorCoordinate,
    sourceOutgoingBoundary,
    displayedStrand: {
      before: displayedStrand(located.strand, currentReverseComplement),
      after: strandAfter
    },
    provenance: {
      schema: 1,
      recordKey: identity.recordKey,
      biologicalFeatureId: identity.biologicalFeatureId,
      placement,
      anchor: placement === 'anchor' ? anchor : null,
      offsetBp: intent.offsetBp,
      orientForward: intent.orientForward
    }
  };
};
