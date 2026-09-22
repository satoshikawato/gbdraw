import assert from 'node:assert/strict';
import test from 'node:test';

import { resolveFeatureAnchor } from '../../gbdraw/web/js/app/record-display/feature-anchor.js';

const identity = { recordKey: 'record-1', biologicalFeatureId: 'feature-1' };

const inputFor = ({
  strand = '+',
  parts = [{ start: 10, end: 16, strand }],
  profile = {},
  intent = {},
  ...record
} = {}) => ({
  recordLength: 100,
  effectiveCircular: true,
  cropped: false,
  currentReverseComplement: false,
  identity,
  parts,
  profile: {
    precision: 'exact',
    operator: parts.length === 1 ? 'single' : 'join',
    partOrder: strand === 'unstranded' ? 'source-forward' : 'biological',
    strand,
    ...profile
  },
  intent: {
    placement: 'anchor',
    anchor: 'five-prime',
    offsetBp: 0,
    orientForward: false,
    ...intent
  },
  ...record
});

const resolve = (options = {}) => resolveFeatureAnchor(inputFor(options));
const coordinate = (options = {}) => {
  const result = resolve(options);
  assert.equal(result.eligibility.enabled, true, result.eligibility.message);
  return result.startCoordinate;
};

test('plus and minus anchors use covered biological traversal', () => {
  assert.equal(coordinate(), 11);
  assert.equal(coordinate({ intent: { anchor: 'midpoint' } }), 13);
  assert.equal(coordinate({ intent: { anchor: 'three-prime' } }), 16);

  assert.equal(coordinate({ strand: '-' }), 16);
  assert.equal(coordinate({ strand: '-', intent: { anchor: 'midpoint' } }), 14);
  assert.equal(coordinate({ strand: '-', intent: { anchor: 'three-prime' } }), 11);
});

test('signed strand-relative offsets wrap in both directions', () => {
  assert.equal(coordinate({ intent: { offsetBp: 95 } }), 6);
  assert.equal(coordinate({ intent: { offsetBp: -12 } }), 99);
  assert.equal(coordinate({ strand: '-', intent: { offsetBp: 20 } }), 96);
  assert.equal(coordinate({ strand: '-', intent: { offsetBp: -90 } }), 6);
});

test('multipart midpoint excludes gaps and uses the earlier covered base for even length', () => {
  assert.equal(coordinate({
    parts: [{ start: 0, end: 3, strand: '+' }, { start: 10, end: 12, strand: '+' }],
    intent: { anchor: 'midpoint' }
  }), 3);
  assert.equal(coordinate({
    parts: [{ start: 0, end: 3, strand: '+' }, { start: 10, end: 13, strand: '+' }],
    intent: { anchor: 'midpoint' }
  }), 3);
  assert.equal(coordinate({
    parts: [{ start: 90, end: 100, strand: '+' }, { start: 0, end: 5, strand: '+' }],
    intent: { anchor: 'midpoint' }
  }), 98);
  assert.equal(coordinate({
    strand: '-',
    parts: [{ start: 0, end: 5, strand: '-' }, { start: 90, end: 100, strand: '-' }],
    intent: { anchor: 'midpoint' }
  }), 98);
});

test('orientation is absolute, preserves current state when off, and is idempotent', () => {
  const preserved = resolve({ currentReverseComplement: true });
  assert.equal(preserved.reverseComplement, true);
  assert.deepEqual(preserved.displayedStrand, { before: '-', after: '-' });

  const forward = resolve({ currentReverseComplement: true, intent: { orientForward: true } });
  assert.equal(forward.reverseComplement, false);
  assert.deepEqual(forward.displayedStrand, { before: '-', after: '+' });
  const repeated = resolveFeatureAnchor({ ...inputFor({ intent: { orientForward: true } }),
    currentReverseComplement: forward.reverseComplement });
  assert.equal(repeated.reverseComplement, forward.reverseComplement);
  assert.equal(repeated.startCoordinate, forward.startCoordinate);
  assert.equal(repeated.displayedStrand.after, forward.displayedStrand.after);
  assert.deepEqual(repeated.provenance, forward.provenance);

  const reverseForward = resolve({ strand: '-', intent: { orientForward: true } });
  assert.equal(reverseForward.reverseComplement, true);
  assert.equal(reverseForward.displayedStrand.after, '+');
});

test('feature-end is the outgoing boundary rather than the 3-prime base', () => {
  const plusThreePrime = coordinate({ intent: { anchor: 'three-prime' } });
  const plusEnd = resolve({ intent: { placement: 'feature-end', anchor: null } });
  assert.equal(plusThreePrime, 16);
  assert.equal(plusEnd.sourceOutgoingBoundary, 17);
  assert.equal(plusEnd.startCoordinate, 17);

  const minusThreePrime = coordinate({ strand: '-', intent: { anchor: 'three-prime', orientForward: true } });
  const minusEnd = resolve({ strand: '-', intent: { placement: 'feature-end', anchor: null, orientForward: true } });
  assert.equal(minusThreePrime, 11);
  assert.equal(minusEnd.sourceOutgoingBoundary, 10);
  assert.equal(minusEnd.startCoordinate, 10);
});

test('unstranded exact paths expose source-forward midpoint and offset wording', () => {
  const result = resolve({
    strand: 'unstranded',
    parts: [{ start: 10, end: 13, strand: 'undefined' }, { start: 20, end: 23, strand: 'undefined' }],
    intent: { anchor: 'midpoint', offsetBp: 2 }
  });
  assert.equal(result.eligibility.enabled, true);
  assert.equal(result.sourceAnchorCoordinate, 13);
  assert.equal(result.startCoordinate, 15);
  assert.equal(result.capabilities.offset.basis, 'source-forward');
  assert.equal(result.capabilities.anchors['five-prime'].code, 'biological-direction-unavailable');
  assert.equal(result.capabilities.orientForward.code, 'orientation-strand-unavailable');
  assert.equal(result.capabilities.featureEnd.enabled, true);
});

test('ambiguous, fuzzy, and unordered locations disable only unsafe operations with reasons', () => {
  const mixed = resolve({
    strand: 'mixed',
    parts: [{ start: 1, end: 3, strand: '+' }, { start: 5, end: 7, strand: '-' }],
    profile: { partOrder: 'ambiguous' },
    intent: { anchor: 'midpoint' }
  });
  assert.equal(mixed.eligibility.code, 'location-traversal-ambiguous');
  assert.equal(mixed.capabilities.orientForward.code, 'orientation-strand-unavailable');

  const fuzzy = resolve({ profile: { precision: 'fuzzy' } });
  assert.equal(fuzzy.eligibility.code, 'location-fuzzy');
  assert.equal(fuzzy.capabilities.orientForward.enabled, true);

  const legacyUnavailable = resolve({ profile: {
    precision: 'unavailable', operator: 'unknown', partOrder: 'ambiguous', strand: '+'
  } });
  assert.equal(legacyUnavailable.eligibility.code, 'feature-metadata-refresh-required');
  assert.match(legacyUnavailable.eligibility.message, /Generate again/);

  const ordered = resolve({ profile: { operator: 'order', partOrder: 'ambiguous' }, intent: { anchor: 'midpoint' } });
  assert.equal(ordered.eligibility.code, 'location-order-ambiguous');
  const unknown = resolve({ profile: { operator: 'unknown', partOrder: 'ambiguous' }, intent: { anchor: 'midpoint' } });
  assert.equal(unknown.eligibility.code, 'location-operator-unknown');

  const repeatedWrap = resolve({
    strand: 'unstranded',
    parts: [
      { start: 70, end: 75, strand: 'undefined' },
      { start: 40, end: 45, strand: 'undefined' },
      { start: 10, end: 15, strand: 'undefined' }
    ],
    intent: { anchor: 'midpoint' }
  });
  assert.equal(repeatedWrap.eligibility.code, 'location-source-order-invalid');
});

test('invalid inputs, zero length, non-circular, and cropped records fail explicitly', () => {
  assert.equal(resolve({ intent: { offsetBp: 1.5 } }).eligibility.code, 'offset-invalid');
  assert.equal(resolve({ intent: { offsetBp: Number.MAX_SAFE_INTEGER + 1 } }).eligibility.code, 'offset-invalid');
  assert.equal(resolve({ recordLength: 0 }).eligibility.code, 'record-length-unavailable');
  assert.equal(resolve({ recordLength: null }).eligibility.code, 'record-length-unavailable');
  assert.equal(resolve({ effectiveCircular: false }).eligibility.code, 'record-not-circular');
  assert.equal(resolve({ cropped: true }).eligibility.code, 'record-cropped');
});

test('resolver does not mutate deeply frozen input and returns deterministic provenance', () => {
  const input = inputFor({ intent: { anchor: 'midpoint', offsetBp: -4, orientForward: true } });
  const deepFreeze = (value) => {
    if (value && typeof value === 'object') {
      Object.values(value).forEach(deepFreeze);
      Object.freeze(value);
    }
    return value;
  };
  deepFreeze(input);

  const first = resolveFeatureAnchor(input);
  const second = resolveFeatureAnchor(input);

  assert.deepEqual(first, second);
  assert.deepEqual(first.provenance, {
    schema: 1,
    recordKey: 'record-1',
    biologicalFeatureId: 'feature-1',
    placement: 'anchor',
    anchor: 'midpoint',
    offsetBp: -4,
    orientForward: true
  });
});
