// Owner decision Q3 = A (design Q4 3.4, 6.3): a successful Generate that
// replaced a source removes only the edits Python reported `unresolved` for a
// replaced record, and the edits of records the request dropped; edits of
// features outside the crop stay dormant, and edits of features the new source
// still has stay, Feature placement rows included.
import assert from 'node:assert/strict';
import test from 'node:test';
import {
  countUnresolvedFeatureEdits,
  pruneUnmatchedFeatureOverrides,
  removeUnresolvedFeatureEdits
} from '../../gbdraw/web/js/services/feature-visibility.js';

// The rows of the Linear drawing (R2; PR-1: each mode's drawing holds its own
// rows, keyed by the identity pair).
const key = (recordKey, featureId) => JSON.stringify([recordKey, featureId]);
const row = (recordKey, biologicalFeatureId, fields) => ({
  recordKey,
  biologicalFeatureId,
  featureVisibility: null,
  labelVisibility: null,
  labelText: null,
  labelSourceText: null,
  ...fields
});
const placement = (recordKey, biologicalFeatureId) => ({
  recordKey, biologicalFeatureId, placement: { kind: 'main' }
});
const notice = (recordKey, biologicalFeatureId, status, kinds) => ({
  recordKey, biologicalFeatureId, status, kinds, resultIndex: 0
});

test('source replacement removes only edits whose feature the new source does not have', () => {
  const featureOverrides = {
    [key('seq-a', 'gone')]: row('seq-a', 'gone', { featureVisibility: 'off', labelText: 'GONE' }),
    [key('seq-a', 'kept')]: row('seq-a', 'kept', { labelVisibility: 'off' }),
    [key('seq-a', 'cropped')]: row('seq-a', 'cropped', { featureVisibility: 'off' }),
    [key('seq-b', 'other')]: row('seq-b', 'other', { featureVisibility: 'on' })
  };
  const featurePlacementOverrides = {
    [key('seq-a', 'gone')]: placement('seq-a', 'gone'),
    [key('seq-a', 'kept')]: placement('seq-a', 'kept')
  };
  const removed = pruneUnmatchedFeatureOverrides({
    featureOverrides,
    featurePlacementOverrides,
    notices: [
      notice('seq-a', 'gone', 'unresolved', ['placement', 'feature_visibility', 'label_text']),
      notice('seq-a', 'cropped', 'crop_excluded', ['feature_visibility'])
    ],
    replacedRecordKeys: ['seq-a'],
    previousRecords: [{ recordKey: 'seq-a' }, { recordKey: 'seq-b' }],
    currentRecords: [{ recordKey: 'seq-a' }, { recordKey: 'seq-b' }],
    biologicalFeatures: [{ scope: 'linear', record_key: 'seq-a', biological_feature_id: 'kept' }]
  });
  assert.equal(removed, 3);
  assert.deepEqual(Object.keys(featureOverrides).sort(), [
    key('seq-a', 'cropped'), key('seq-a', 'kept'), key('seq-b', 'other')
  ].sort());
  assert.deepEqual(Object.keys(featurePlacementOverrides), [key('seq-a', 'kept')]);
});

test('a record the request dropped loses its edits; the other drawing keeps its own (R2)', () => {
  const featureOverrides = {
    [key('seq-a', 'f1')]: row('seq-a', 'f1', { featureVisibility: 'off', labelText: 'A' })
  };
  // The Circular drawing's rows on the same record key are its own map.
  const circularOverrides = {
    [key('seq-a', 'f1')]: row('seq-a', 'f1', { labelVisibility: 'on' }),
    [key('circular-x', 'f1')]: row('circular-x', 'f1', { labelVisibility: 'on' })
  };
  const circularBefore = structuredClone(circularOverrides);
  const featurePlacementOverrides = { [key('seq-a', 'f1')]: placement('seq-a', 'f1') };
  const removed = pruneUnmatchedFeatureOverrides({
    featureOverrides,
    featurePlacementOverrides,
    notices: [],
    replacedRecordKeys: ['seq-a'],
    previousRecords: [{ recordKey: 'seq-a' }, { recordKey: 'seq-b' }],
    currentRecords: [{ recordKey: 'seq-b' }],
    biologicalFeatures: []
  });
  assert.equal(removed, 3);
  assert.deepEqual(featureOverrides, {});
  assert.deepEqual(featurePlacementOverrides, {});
  assert.deepEqual(circularOverrides, circularBefore);
});

test('a kept label source text goes with a feature the replaced source lost', () => {
  const featureOverrides = {
    [key('seq-a', 'gone')]: row('seq-a', 'gone', { labelSourceText: 'alpha' }),
    [key('seq-a', 'kept')]: row('seq-a', 'kept', { labelSourceText: 'beta' })
  };
  assert.equal(pruneUnmatchedFeatureOverrides({
    featureOverrides,
    replacedRecordKeys: ['seq-a'],
    biologicalFeatures: [{ scope: 'linear', record_key: 'seq-a', biological_feature_id: 'kept' }]
  }), 0);
  assert.deepEqual(Object.keys(featureOverrides), [key('seq-a', 'kept')]);
});

test('without a replaced source unresolved edits stay until the user removes them', () => {
  const featureOverrides = { [key('seq-a', 'gone')]: row('seq-a', 'gone', { labelText: 'X', featureVisibility: 'off' }) };
  const featurePlacementOverrides = { [key('seq-a', 'gone')]: placement('seq-a', 'gone') };
  const notices = [notice('seq-a', 'gone', 'unresolved', ['placement', 'feature_visibility', 'label_text'])];
  assert.equal(pruneUnmatchedFeatureOverrides({
    featureOverrides, featurePlacementOverrides, notices, replacedRecordKeys: []
  }), 0);
  assert.equal(countUnresolvedFeatureEdits({ featureOverrides, featurePlacementOverrides, notices }), 3);
  assert.equal(removeUnresolvedFeatureEdits({ featureOverrides, featurePlacementOverrides, notices }), 3);
  assert.deepEqual([featureOverrides, featurePlacementOverrides], [{}, {}]);
  assert.equal(countUnresolvedFeatureEdits({ featureOverrides, featurePlacementOverrides, notices }), 0);
});

// OV-84 (Owner-delegated 2026-10-07): a source replacement retires the
// per-feature strokes of features the new source does not have, as it does
// their label edits; a stroke of a feature the new source still has, and a
// stroke on a record the Generate did not replace, stay.
test('source replacement retires the strokes of features the new source lost (OV-84)', () => {
  const stroke = (recordKey, featureId) => `${recordKey}\0${featureId}`;
  const featureStrokeOverrides = {
    [stroke('seq-a', 'gone')]: { strokeColor: '#111111', strokeWidth: 2 },
    [stroke('seq-a', 'kept')]: { strokeColor: '#222222', strokeWidth: 3 },
    [stroke('seq-b', 'other')]: { strokeColor: '#333333', strokeWidth: 1 }
  };
  const removed = pruneUnmatchedFeatureOverrides({
    featureStrokeOverrides,
    replacedRecordKeys: ['seq-a'],
    previousRecords: [{ recordKey: 'seq-a' }, { recordKey: 'seq-b' }],
    currentRecords: [{ recordKey: 'seq-a' }, { recordKey: 'seq-b' }],
    biologicalFeatures: [{ scope: 'linear', record_key: 'seq-a', biological_feature_id: 'kept' }]
  });
  assert.equal(removed, 1);
  assert.deepEqual(Object.keys(featureStrokeOverrides).sort(), [stroke('seq-a', 'kept'), stroke('seq-b', 'other')].sort());
  // Without a replaced source, every stroke stays dormant until its feature returns.
  assert.equal(pruneUnmatchedFeatureOverrides({ featureStrokeOverrides, replacedRecordKeys: [] }), 0);
  assert.equal(Object.keys(featureStrokeOverrides).length, 2);
});
