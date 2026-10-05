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
} from '../../gbdraw/web/js/app/feature-visibility.js';

// The rows of a Linear request; R2 scopes each draft row to its mode.
const key = (recordKey, featureId, scope = 'linear') => JSON.stringify([scope, recordKey, featureId]);
const row = (recordKey, biologicalFeatureId, fields) => ({
  scope: 'linear',
  recordKey,
  biologicalFeatureId,
  featureVisibility: null,
  labelVisibility: null,
  labelText: null,
  labelSourceText: null,
  ...fields
});
const placement = (recordKey, biologicalFeatureId) => ({
  scope: 'linear', recordKey, biologicalFeatureId, placement: { kind: 'main' }
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
    scope: 'linear',
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

test('a record the request dropped loses its edits; another mode keeps its own (R2)', () => {
  const featureOverrides = {
    [key('seq-a', 'f1')]: row('seq-a', 'f1', { featureVisibility: 'off', labelText: 'A' }),
    [key('seq-a', 'f1', 'circular')]: row('seq-a', 'f1', { scope: 'circular', labelVisibility: 'on' }),
    [key('circular-x', 'f1', 'circular')]: row('circular-x', 'f1', { scope: 'circular', labelVisibility: 'on' })
  };
  const featurePlacementOverrides = { [key('seq-a', 'f1')]: placement('seq-a', 'f1') };
  const removed = pruneUnmatchedFeatureOverrides({
    featureOverrides,
    featurePlacementOverrides,
    notices: [],
    scope: 'linear',
    replacedRecordKeys: ['seq-a'],
    previousRecords: [{ recordKey: 'seq-a' }, { recordKey: 'seq-b' }],
    currentRecords: [{ recordKey: 'seq-b' }],
    biologicalFeatures: []
  });
  assert.equal(removed, 3);
  assert.deepEqual(Object.keys(featureOverrides), [key('seq-a', 'f1', 'circular'), key('circular-x', 'f1', 'circular')]);
  assert.deepEqual(featurePlacementOverrides, {});
});

test('a kept label source text goes with a feature the replaced source lost', () => {
  const featureOverrides = {
    [key('seq-a', 'gone')]: row('seq-a', 'gone', { labelSourceText: 'alpha' }),
    [key('seq-a', 'kept')]: row('seq-a', 'kept', { labelSourceText: 'beta' })
  };
  assert.equal(pruneUnmatchedFeatureOverrides({
    featureOverrides,
    scope: 'linear',
    replacedRecordKeys: ['seq-a'],
    biologicalFeatures: [{ scope: 'linear', record_key: 'seq-a', biological_feature_id: 'kept' }]
  }), 0);
  assert.deepEqual(Object.keys(featureOverrides), [key('seq-a', 'kept')]);
});

test('without a replaced source unresolved edits stay until the user removes them', () => {
  const featureOverrides = { [key('seq-a', 'gone')]: row('seq-a', 'gone', { labelText: 'X', featureVisibility: 'off' }) };
  const featurePlacementOverrides = { [key('seq-a', 'gone')]: placement('seq-a', 'gone') };
  const notices = [notice('seq-a', 'gone', 'unresolved', ['placement', 'feature_visibility', 'label_text'])];
  const scope = 'linear';
  assert.equal(pruneUnmatchedFeatureOverrides({
    featureOverrides, featurePlacementOverrides, notices, scope, replacedRecordKeys: []
  }), 0);
  assert.equal(countUnresolvedFeatureEdits({ featureOverrides, featurePlacementOverrides, notices, scope }), 3);
  assert.equal(removeUnresolvedFeatureEdits({ featureOverrides, featurePlacementOverrides, notices, scope }), 3);
  assert.deepEqual([featureOverrides, featurePlacementOverrides], [{}, {}]);
  assert.equal(countUnresolvedFeatureEdits({ featureOverrides, featurePlacementOverrides, notices, scope }), 0);
});
