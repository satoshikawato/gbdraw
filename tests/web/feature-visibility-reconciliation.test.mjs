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
import { FEATURE_EDITS_SLICE } from '../../gbdraw/web/js/services/drawing/feature-edits.js';
import { PAINT_SLICE } from '../../gbdraw/web/js/services/drawing/paint.js';
import { ANNOTATION_SETS_SLICE } from '../../gbdraw/web/js/services/drawing/annotation-sets.js';

// The rows of the Linear drawing (R2; PR-1: each mode's drawing holds its own
// rows, keyed by the identity pair).
// A drawing's feature-edits slice over its draft maps, so a writer writes them.
const editsOf = (featureOverrides, featurePlacementOverrides) => (
  FEATURE_EDITS_SLICE.create({ featureOverrides, featurePlacementOverrides }).handle
);
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
    featureEdits: editsOf(featureOverrides, featurePlacementOverrides),
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
    featureEdits: editsOf(featureOverrides, featurePlacementOverrides),
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
    featureEdits: editsOf(featureOverrides),
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
    featureEdits: editsOf(featureOverrides, featurePlacementOverrides), notices, replacedRecordKeys: []
  }), 0);
  assert.equal(countUnresolvedFeatureEdits({ featureEdits: editsOf(featureOverrides, featurePlacementOverrides), notices }), 3);
  assert.equal(removeUnresolvedFeatureEdits({ featureEdits: editsOf(featureOverrides, featurePlacementOverrides), notices }), 3);
  assert.deepEqual([featureOverrides, featurePlacementOverrides], [{}, {}]);
  assert.equal(countUnresolvedFeatureEdits({ featureEdits: editsOf(featureOverrides, featurePlacementOverrides), notices }), 0);
});

// OV-84 (Owner-delegated 2026-10-07): a source replacement retires the
// per-feature strokes of features the new source does not have, as it does
// their label edits; a stroke of a feature the new source still has, and a
// stroke on a record the Generate did not replace, stay.
test('source replacement retires the strokes of features the new source lost (OV-84)', () => {
  const stroke = (recordKey, featureId) => `${recordKey}\0${featureId}`;
  const paint = PAINT_SLICE.create({ featureStrokeOverrides: {
    [stroke('seq-a', 'gone')]: { strokeColor: '#111111', strokeWidth: 2 },
    [stroke('seq-a', 'kept')]: { strokeColor: '#222222', strokeWidth: 3 },
    [stroke('seq-b', 'other')]: { strokeColor: '#333333', strokeWidth: 1 }
  } });
  const { featureStrokeOverrides } = paint.view;
  const removed = pruneUnmatchedFeatureOverrides({
    featureEdits: editsOf(), paint: paint.handle,
    replacedRecordKeys: ['seq-a'],
    previousRecords: [{ recordKey: 'seq-a' }, { recordKey: 'seq-b' }],
    currentRecords: [{ recordKey: 'seq-a' }, { recordKey: 'seq-b' }],
    biologicalFeatures: [{ scope: 'linear', record_key: 'seq-a', biological_feature_id: 'kept' }]
  });
  assert.equal(removed, 1);
  assert.deepEqual(Object.keys(featureStrokeOverrides).sort(), [stroke('seq-a', 'kept'), stroke('seq-b', 'other')].sort());
  // Without a replaced source, every stroke stays dormant until its feature returns.
  assert.equal(pruneUnmatchedFeatureOverrides({ featureEdits: editsOf(), paint: paint.handle, replacedRecordKeys: [] }), 0);
  assert.equal(Object.keys(featureStrokeOverrides).length, 2);
});

test('an OFF record is not dropped: its edits wait for it (record selection D-05)', () => {
  const featureOverrides = {
    [key('seq-a', 'f1')]: row('seq-a', 'f1', { featureVisibility: 'off' }),
    [key('seq-b', 'f2')]: row('seq-b', 'f2', { labelVisibility: 'off' })
  };
  const featurePlacementOverrides = { [key('seq-a', 'f1')]: placement('seq-a', 'f1') };
  const featureStrokeOverrides = { 'seq-a\0stable-1': { width: 2 } };
  const before = structuredClone({ featureOverrides, featurePlacementOverrides, featureStrokeOverrides });
  const edits = FEATURE_EDITS_SLICE.create({ featureOverrides, featurePlacementOverrides });
  const paint = PAINT_SLICE.create({ featureStrokeOverrides });
  const removed = pruneUnmatchedFeatureOverrides({
    featureEdits: edits.handle,
    paint: paint.handle,
    notices: [],
    replacedRecordKeys: ['seq-a', 'seq-b'],
    previousRecords: [{ recordKey: 'seq-a' }, { recordKey: 'seq-b' }],
    currentRecords: [{ recordKey: 'seq-b' }],
    retainedRecords: [{ recordKey: 'seq-a', cardinality: 'all' }],
    biologicalFeatures: [{ scope: 'linear', record_key: 'seq-b', biological_feature_id: 'f2' }]
  });
  assert.equal(removed, 0);
  assert.deepEqual({
    featureOverrides: edits.view.featureOverrides,
    featurePlacementOverrides: edits.view.featurePlacementOverrides,
    featureStrokeOverrides: paint.view.featureStrokeOverrides
  }, before);
});

// OV-1949 (record-identity-20261027 DESIGN A5): every successful Generate
// moves the rows of each record whose key changed; the removal of a
// source-replacing Generate reads the keys before the move.
const movesOf = (entries) => new Map(entries);

test('a Generate moves a record\'s rows, strokes and annotation targets to its new key', () => {
  const featureOverrides = {
    [key('record-2', 'f1')]: row('record-2', 'f1', { featureVisibility: 'off', labelText: 'moved' }),
    [key('record-1', 'f1')]: row('record-1', 'f1', { labelText: 'stays' })
  };
  const featurePlacementOverrides = { [key('record-2', 'f1')]: placement('record-2', 'f1') };
  const featureStrokeOverrides = { 'record-2\0f1': { width: 2 }, 'record-1\0f1': { width: 3 } };
  const edits = FEATURE_EDITS_SLICE.create({ featureOverrides, featurePlacementOverrides });
  const paint = PAINT_SLICE.create({ featureStrokeOverrides });
  const annotations = ANNOTATION_SETS_SLICE.create({ annotationSets: [{ id: 'set', annotations: [
    { target: { kind: 'featureIdentity', recordKey: 'record-2', biologicalFeatureId: 'f1' } },
    { target: { kind: 'featureIdentity', recordKey: 'record-1', biologicalFeatureId: 'f1' } }
  ] }] });
  const removed = pruneUnmatchedFeatureOverrides({
    featureEdits: edits.handle,
    paint: paint.handle,
    annotationSets: annotations.handle,
    annotationSetRows: annotations.view.annotationSets,
    currentRecords: [{ recordKey: 'circular-B-#2' }],
    moves: movesOf([['record-2', 'circular-B-#2']])
  });
  assert.equal(removed, 0);
  assert.deepEqual(edits.view.featureOverrides, {
    [key('circular-B-#2', 'f1')]: row('circular-B-#2', 'f1', { featureVisibility: 'off', labelText: 'moved' }),
    [key('record-1', 'f1')]: row('record-1', 'f1', { labelText: 'stays' })
  });
  assert.deepEqual(edits.view.featurePlacementOverrides, { [key('circular-B-#2', 'f1')]: placement('circular-B-#2', 'f1') });
  assert.deepEqual(paint.view.featureStrokeOverrides, { 'circular-B-#2\0f1': { width: 2 }, 'record-1\0f1': { width: 3 } });
  assert.deepEqual(annotations.view.annotationSets[0].annotations.map((item) => item.target.recordKey), ['circular-B-#2', 'record-1']);
});

test('rows meeting under one key: the record\'s own row wins field by field over an alias row', () => {
  const featureOverrides = {
    [key('circular-B-#2', 'f1')]: row('circular-B-#2', 'f1', { labelText: 'own', featureVisibility: null }),
    [key('record-2', 'f1')]: row('record-2', 'f1', { labelText: 'alias', featureVisibility: 'off' })
  };
  const edits = FEATURE_EDITS_SLICE.create({ featureOverrides, featurePlacementOverrides: {} });
  // The record's own key moves first; its alias second (`recordKeyMoves` order).
  pruneUnmatchedFeatureOverrides({
    featureEdits: edits.handle,
    moves: movesOf([['circular-B-#2', 'record-2-canvas'], ['record-2', 'record-2-canvas']])
  });
  assert.deepEqual(edits.view.featureOverrides, {
    [key('record-2-canvas', 'f1')]: row('record-2-canvas', 'f1', { labelText: 'own', featureVisibility: 'off' })
  });
});

test('a permutation moves every row together from one snapshot', () => {
  const featureOverrides = {
    [key('record-1', 'f1')]: row('record-1', 'f1', { labelText: 'R3' }),
    [key('record-2', 'f1')]: row('record-2', 'f1', { labelText: 'R1' })
  };
  const edits = FEATURE_EDITS_SLICE.create({ featureOverrides, featurePlacementOverrides: {} });
  pruneUnmatchedFeatureOverrides({ featureEdits: edits.handle, moves: movesOf([['record-1', 'record-3'], ['record-2', 'record-1']]) });
  assert.deepEqual(Object.values(edits.view.featureOverrides).map((item) => [item.recordKey, item.labelText]).sort(),
    [['record-1', 'R1'], ['record-3', 'R3']]);
});

test('A5: a source-replacing Generate drops a record the request no longer draws, before the move', () => {
  // Card U drew every record (U:1, U:2) and now draws record 1 (key U); U's source is replaced.
  const featureOverrides = {
    [key('U:1', 'f1')]: row('U:1', 'f1', { labelText: 'one' }),
    [key('U:2', 'f2')]: row('U:2', 'f2', { labelText: 'two' })
  };
  const replacing = FEATURE_EDITS_SLICE.create({ featureOverrides: structuredClone(featureOverrides), featurePlacementOverrides: {} });
  pruneUnmatchedFeatureOverrides({
    featureEdits: replacing.handle,
    replacedRecordKeys: ['U:1', 'U:2'],
    previousRecords: [{ recordKey: 'U', cardinality: 'all' }],
    currentRecords: [{ recordKey: 'U', cardinality: 'exactly_one' }],
    // The new catalog names record 1's feature under the new key: it is not gone.
    biologicalFeatures: [{ record_key: 'U', biological_feature_id: 'f1' }],
    moves: movesOf([['U:1', 'U']])
  });
  assert.deepEqual(replacing.view.featureOverrides, { [key('U', 'f1')]: row('U', 'f1', { labelText: 'one' }) });
  // Without a source replacement, record 2's rows wait under its identity.
  const keeping = FEATURE_EDITS_SLICE.create({ featureOverrides: structuredClone(featureOverrides), featurePlacementOverrides: {} });
  pruneUnmatchedFeatureOverrides({
    featureEdits: keeping.handle,
    currentRecords: [{ recordKey: 'U', cardinality: 'exactly_one' }],
    moves: movesOf([['U:1', 'U']])
  });
  assert.deepEqual(keeping.view.featureOverrides, {
    [key('U', 'f1')]: row('U', 'f1', { labelText: 'one' }),
    [key('U:2', 'f2')]: row('U:2', 'f2', { labelText: 'two' })
  });
});

test('M2: an OFF card keeps its rows, under its uid or its identity, when another card\'s source is replaced', () => {
  const featureOverrides = {
    [key('V', 'f1')]: row('V', 'f1', { labelText: 'drawn alone' }),
    [key('V:2', 'f2')]: row('V:2', 'f2', { labelText: 'waiting' }),
    [key('U', 'f3')]: row('U', 'f3', { labelText: 'u' })
  };
  const edits = FEATURE_EDITS_SLICE.create({ featureOverrides, featurePlacementOverrides: {} });
  const removed = pruneUnmatchedFeatureOverrides({
    featureEdits: edits.handle,
    replacedRecordKeys: ['U'],
    previousRecords: [{ recordKey: 'U', cardinality: 'all' }, { recordKey: 'V', cardinality: 'exactly_one' }, { recordKey: 'V:2', cardinality: 'exactly_one' }],
    currentRecords: [{ recordKey: 'U', cardinality: 'all' }],
    retainedRecords: [{ recordKey: 'V', cardinality: 'all' }],
    biologicalFeatures: [{ record_key: 'U', biological_feature_id: 'f3' }],
    // V went OFF: it drew record 2 alone, so its rows wait under V:2.
    moves: movesOf([['V', 'V:2']])
  });
  assert.equal(removed, 0);
  assert.deepEqual(Object.keys(edits.view.featureOverrides).sort(), [key('U', 'f3'), key('V:2', 'f1'), key('V:2', 'f2')].sort());
});
