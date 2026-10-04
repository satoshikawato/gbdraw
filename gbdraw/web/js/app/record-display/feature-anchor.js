const available = (extra = {}) => ({ enabled: true, code: null, message: '', ...extra });
const unavailable = (code, message, extra = {}) => ({ enabled: false, code, message, ...extra });

const firstUnavailable = (...operations) => operations.find((operation) => !operation.enabled) || null;

const wrapCoordinate = (coordinate, length) => 1 + ((((coordinate - 1) % length) + length) % length);

const displayedStrand = (sourceStrand, reverseComplement) => {
  if (!['+', '-'].includes(sourceStrand)) return sourceStrand;
  if (!reverseComplement) return sourceStrand;
  return sourceStrand === '+' ? '-' : '+';
};

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

const partStrandMatches = (partStrand, profileStrand) => {
  if (profileStrand === 'mixed') return true;
  if (profileStrand === 'unstranded') {
    return !['+', '-'].includes(partStrand);
  }
  return partStrand === profileStrand;
};

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
  const sourceOrder = [...parts].sort((left, right) => left.start - right.start || left.end - right.end);
  if (sourceOrder.some((part, index) => index > 0 && part.start < sourceOrder[index - 1].end)) {
    return unavailable('location-parts-overlap', 'Overlapping feature parts cannot be traversed safely.');
  }
  if (profile.partOrder === 'source-forward') {
    const wraps = parts.reduce((count, part, index) => count
      + (index > 0 && part.start < parts[index - 1].start ? 1 : 0), 0);
    if (wraps > 1) {
      return unavailable('location-source-order-invalid', 'Feature parts do not form one source-forward circular path.');
    }
  }
  return available();
};

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

const intentCapability = (intent) => {
  if (!intent || !['anchor', 'feature-end'].includes(intent.placement)) {
    return unavailable('placement-invalid', 'Choose an anchor or feature-end placement.');
  }
  if (intent.placement === 'anchor'
    && !['five-prime', 'midpoint', 'three-prime'].includes(intent.anchor)) {
    return unavailable('anchor-invalid', 'Choose a supported feature anchor.');
  }
  if (intent.placement === 'feature-end' && intent.anchor !== null) {
    return unavailable('anchor-not-applicable', 'Feature-end placement does not use an anchor.');
  }
  if (!Number.isSafeInteger(intent.offsetBp)) {
    return unavailable('offset-invalid', 'Offset must be a safe integer number of base pairs.');
  }
  if (intent.orientForward !== true && intent.orientForward !== false) {
    return unavailable('orient-forward-invalid', 'Orient-forward intent must be true or false.');
  }
  return available();
};

const identityCapability = (identity) => {
  if (!identity || typeof identity.recordKey !== 'string' || !identity.recordKey
    || identity.recordKey.includes('\0') || typeof identity.biologicalFeatureId !== 'string'
    || !identity.biologicalFeatureId || identity.biologicalFeatureId.includes('\0')) {
    return unavailable('feature-identity-unavailable', 'Feature rotation requires a stable record and feature identity.');
  }
  return available();
};

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
  const partsState = profileState.enabled
    ? partsCapability(parts, profile, recordLength) : profileState;
  const positional = positionalCapability(record, profileState, partsState, profile || {});

  const midpoint = positional;
  const biologicalEnds = !positional.enabled ? positional
    : !['+', '-'].includes(profile.strand)
      ? unavailable('biological-direction-unavailable', 'Biological 5-prime and 3-prime ends require a known strand.')
      : available();
  const offsetBasis = profile?.strand === 'unstranded' ? 'source-forward' : 'strand';
  const offset = !positional.enabled ? { ...positional, basis: offsetBasis }
    : available({ basis: offsetBasis });
  const orientForward = !record.enabled ? record
    : !profileState.enabled ? profileState
      : !['+', '-'].includes(profile.strand)
        ? unavailable('orientation-strand-unavailable', 'Orient forward requires one known feature strand.')
        : available();
  const featureEnd = !positional.enabled ? positional
    : profile.strand === 'unstranded' && profile.partOrder !== 'source-forward'
      ? unavailable('feature-end-direction-unavailable', 'Feature-end placement requires a known display traversal path.')
      : available();
  const capabilities = {
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
    : firstUnavailable(capabilities.anchors[intent.anchor], offset);
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

  const reverseComplement = intent.orientForward ? profile.strand === '-' : currentReverseComplement;
  const direction = profile.strand === '-' ? -1 : 1;
  const coveredLength = parts.reduce((total, part) => total + part.end - part.start, 0);
  let sourceAnchorCoordinate = null;
  let sourceOutgoingBoundary = null;
  let baseCoordinate;
  if (intent.placement === 'anchor') {
    const coveredOffset = intent.anchor === 'five-prime' ? 0
      : intent.anchor === 'midpoint' ? Math.floor((coveredLength - 1) / 2)
        : coveredLength - 1;
    sourceAnchorCoordinate = baseAtCoveredOffset(parts, coveredOffset, direction);
    baseCoordinate = sourceAnchorCoordinate;
  } else {
    const displayStep = reverseComplement ? -1 : 1;
    const followsPath = displayStep === direction;
    const lastOffset = followsPath ? coveredLength - 1 : 0;
    const lastIncludedBase = baseAtCoveredOffset(parts, lastOffset, direction);
    sourceOutgoingBoundary = wrapCoordinate(lastIncludedBase + displayStep, recordLength);
    baseCoordinate = sourceOutgoingBoundary;
  }
  const offsetDirection = profile.strand === '-' ? -1 : 1;
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
      before: displayedStrand(profile.strand, currentReverseComplement),
      after: displayedStrand(profile.strand, reverseComplement)
    },
    provenance: {
      schema: 1,
      recordKey: identity.recordKey,
      biologicalFeatureId: identity.biologicalFeatureId,
      placement: intent.placement,
      anchor: intent.placement === 'anchor' ? intent.anchor : null,
      offsetBp: intent.offsetBp,
      orientForward: intent.orientForward
    }
  };
};
